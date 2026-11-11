#20260915 @Yan Zhong
# -*- coding: utf-8 -*-
"""
Multi.Sim v1.0 - Multi-Hazard Simulation | ArcGIS Python Script Tool
Parameters:
  0  = DEM
  1  = RIA prone raster
  2  = GLOF lakes raster
  3  = Landslide prone raster
  4  = Output folder
  5  = tan alpha RIA          (optional, default 0.10)
  6  = tan alpha GLOF         (optional, default 0.05)
  7  = tan alpha LS           (optional, default 0.19)
  8  = BFS max neighbours     (optional, default 2)
  9  = Min patch px RIA       (optional, default 22)
  10 = Min patch px LS        (optional, default 4)
  11 = Infrastructure layers  (optional, multi-value; point / polyline / polygon)
Outputs:
  RIA_HI.tif, RIA_HI_raw.tif, GLOF_HI.tif, GLOF_HI_raw.tif,
  LS_HI.tif, LS_HI_raw.tif, Multi_HI.tif
  <LayerName>_exposure.shp per infrastructure layer
Exposure fields per active hazard (+ Multi):
  EI_*  Float, 0-1 inside the hazard zone, -1 outside
  CL_*  Short, exposure level 0-5
  Polylines also get EI_*_mx (maximum EI along the line); CL_* uses the mean.
  Classes: 0 Negligible | 1 Very low <=0.2 | 2 Low <=0.4
           3 Moderate <=0.6 | 4 High <=0.8 | 5 Very high >0.8
"""
import arcpy
from arcpy.sa import *
import os, sys, traceback
import numpy as np
import math
from scipy.ndimage import binary_dilation, generate_binary_structure, label
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True
DEFAULT_TANA_RIA      = 0.10
DEFAULT_TANA_GLOF     = 0.05
DEFAULT_TANA_LS       = 0.19
DEFAULT_BFS_MAX_N     = 2
DEFAULT_MIN_PATCH_RIA = 22
DEFAULT_MIN_PATCH_LS  = 4
GLOF_DILATION         = 5
FOCAL_ITERATIONS      = 1
# The running reach angle (elev0 - elev)/path_len is not monotonic: where a path
# crosses standing water the elevation stops falling while path_len keeps rising,
# so the angle dips and then recovers. Stopping at the first violation arrests the
# path inside that dip. Instead the walk continues for this many metres below the
# angle and the runout limit is the farthest node that satisfied it.
#   None -> never give up (GLOF: a confined flood does not stop on a lake surface)
#   0.0  -> first-violation rule, i.e. v1 behaviour
RECOVERY_WINDOW_M = {"GLOF": None, "RIA": 3000.0, "LS": 3000.0}
MAX_PATH_M = 300000.0
ALL_NEIGHBOURS = [
    (-1, -1, math.sqrt(2)), (-1, 0, 1.0), (-1, 1, math.sqrt(2)),
    ( 0, -1, 1.0),                         ( 0, 1, 1.0),
    ( 1, -1, math.sqrt(2)), ( 1, 0, 1.0), ( 1, 1, math.sqrt(2)),
]
def get_raster_center(path):
    ext = arcpy.Describe(path).extent
    return (ext.XMin + ext.XMax) / 2.0, (ext.YMin + ext.YMax) / 2.0
def needs_projection(dataset_path):
    sr = arcpy.Describe(dataset_path).spatialReference
    if sr is None or sr.name in ("", "Unknown"):
        return True, "No spatial reference"
    if sr.type == "Geographic":
        return True, "Geographic CRS: {}".format(sr.name)
    unit = sr.linearUnitName.lower()
    if "meter" in unit or "metre" in unit:
        return False, "Projected CRS (meter): {}".format(sr.name)
    return True, "Projected CRS, non-meter ({}): {}".format(sr.linearUnitName, sr.name)
def build_utm_sr(lon, lat):
    zone = int((lon + 180) / 6) + 1
    hemi = "N" if lat >= 0 else "S"
    epsg = (32600 + zone) if lat >= 0 else (32700 + zone)
    for attempt in [lambda: arcpy.SpatialReference(epsg),
                    lambda: arcpy.SpatialReference(
                        "WGS 1984 UTM Zone {}{}".format(zone, hemi))]:
        try:
            sr = attempt(); _ = sr.name; return sr, epsg
        except Exception:
            pass
    sr = arcpy.SpatialReference(); sr.factoryCode = epsg
    return sr, epsg
def align_raster(input_path, dem_work, out_path):
    dem_sr = arcpy.Describe(dem_work).spatialReference
    src_sr = arcpy.Describe(input_path).spatialReference
    if src_sr.factoryCode != dem_sr.factoryCode:
        arcpy.AddMessage("  Reprojecting {} to DEM CRS...".format(
            os.path.basename(input_path)))
        arcpy.ProjectRaster_management(input_path, out_path, dem_sr, "NEAREST")
        return out_path
    return input_path
def flowdir_to_offset(code, cellsize):
    mapping = {
        1:   ( 0,  1, cellsize),
        2:   ( 1,  1, math.sqrt(2) * cellsize),
        4:   ( 1,  0, cellsize),
        8:   ( 1, -1, math.sqrt(2) * cellsize),
        16:  ( 0, -1, cellsize),
        32:  (-1, -1, math.sqrt(2) * cellsize),
        64:  (-1,  0, cellsize),
        128: (-1,  1, math.sqrt(2) * cellsize),
    }
    return mapping.get(code)
def normalize(arr, percentile=99):
    arr       = arr.copy()
    positives = arr[arr > 0]
    if len(positives) == 0:
        return arr.astype(np.float32)
    arr    = np.clip(arr, 0, np.percentile(positives, percentile))
    mn, mx = arr.min(), arr.max()
    if mx > mn:
        return ((arr - mn) / (mx - mn)).astype(np.float32)
    return arr.astype(np.float32)
def save_raster(arr, ll, pixel_size, sr, path):
    r = arcpy.NumPyArrayToRaster(arr, ll, pixel_size, pixel_size, value_to_nodata=0)
    arcpy.DefineProjection_management(r, sr)
    r.save(path)
def get_downslope_neighbours(r, c, dem_arr, cellsize, nrows, ncols, bfs_max_n):
    elev = dem_arr[r, c]
    candidates = []
    for dr, dc, df in ALL_NEIGHBOURS:
        rn, cn = r + dr, c + dc
        if 0 <= rn < nrows and 0 <= cn < ncols:
            elev_n = dem_arr[rn, cn]
            if np.isfinite(elev_n) and elev_n <= elev:
                candidates.append(
                    (rn, cn, df * cellsize, (elev - elev_n) / (df * cellsize)))
    if not candidates:
        return []
    candidates.sort(key=lambda x: x[3], reverse=True)
    return [(rn, cn, dist) for rn, cn, dist, _ in candidates[:bfs_max_n]]
def remove_small_patches(arr, min_patch_px, name):
    out    = arr.copy()
    struct = generate_binary_structure(2, 2)
    lbl, n = label((out > 0).astype(np.uint8), structure=struct)
    removed = 0
    for lid in range(1, n + 1):
        if np.sum(lbl == lid) < min_patch_px:
            out[lbl == lid] = 0
            removed += 1
    arcpy.AddMessage("    {}: {} small patches removed (< {} px).".format(
        name, removed, min_patch_px))
    return out
def prepare_ria_aoi(prone_arr, min_patch_ria):
    return remove_small_patches(prone_arr, min_patch_ria, "RIA AOI")
def prepare_ls_aoi(ls_arr, ria_arr, min_patch_ls):
    arr = remove_small_patches(ls_arr, min_patch_ls, "LS AOI")
    if ria_arr is not None:
        overlap = (arr > 0) & (ria_arr > 0)
        arr[overlap] = 0
        arcpy.AddMessage("    LS AOI: {} pixels removed (overlap with RIA AOI).".format(
            int(np.sum(overlap))))
    return arr
def sim_path(r0, c0, elev0, tana, dem_arr, flow_arr,
             nrows, ncols, pixel_size, bfs_max_n, weight=1.0,
             recovery_window=0.0):
    # Step 1: D8 main path. Every visited node is recorded; the path is truncated
    # afterwards to the farthest node satisfying the reach angle.
    d8_nodes = []
    d8_set   = set()
    r_cur, c_cur = r0, c0
    path_len = 0.0
    last_ok_idx = 0
    l_last_ok   = 0.0
    while True:
        if (r_cur, c_cur) in d8_set:
            break
        d8_set.add((r_cur, c_cur))
        elev_cur = dem_arr[r_cur, c_cur]
        if not np.isfinite(elev_cur):
            break
        d8_nodes.append((r_cur, c_cur, path_len))
        if path_len > 0:
            if (elev0 - elev_cur) / path_len >= tana:
                last_ok_idx = len(d8_nodes) - 1
                l_last_ok   = path_len
            elif recovery_window is not None and \
                    (path_len - l_last_ok) > recovery_window:
                break
        if path_len > MAX_PATH_M:
            break
        offset = flowdir_to_offset(flow_arr[r_cur, c_cur], pixel_size)
        if offset is None:
            break
        dr, dc, dist = offset
        r_next, c_next = r_cur + dr, c_cur + dc
        if not (0 <= r_next < nrows and 0 <= c_next < ncols):
            break
        path_len += dist
        r_cur, c_cur = r_next, c_next
    d8_nodes = d8_nodes[:last_ok_idx + 1]
    if not d8_nodes:
        return {}
    total_len = d8_nodes[-1][2]
    if total_len <= 0.0:
        total_len = pixel_size
    # Step 2: BFS lateral extent. Each branch carries its cumulative distance from
    # its originating D8 node and may not exceed total_len, so lateral spread
    # cannot become an unconstrained downslope flood-fill.
    bfs_extent     = set()
    spread_visited = {(rd, cd) for (rd, cd, _pl) in d8_nodes}
    for (rd, cd, pl) in d8_nodes:
        queue = [(rd, cd, pl)]
        qi = 0
        while qi < len(queue):
            rs, cs, cur_dist = queue[qi]; qi += 1
            bfs_extent.add((rs, cs))
            for rn, cn, step_dist in get_downslope_neighbours(rs, cs, dem_arr,
                                                       pixel_size, nrows, ncols,
                                                       bfs_max_n):
                new_dist = cur_dist + step_dist
                if new_dist > total_len:
                    continue
                if (rn, cn) not in spread_visited:
                    spread_visited.add((rn, cn))
                    queue.append((rn, cn, new_dist))
    # Step 3: nearest D8 node, quadratic decay
    d8_ei   = {(rd, cd): weight * (1.0 - min(pl / total_len, 1.0)) ** 2
               for (rd, cd, pl) in d8_nodes}
    d8_list = [(rd, cd) for (rd, cd, _pl) in d8_nodes]
    d8_arr  = np.array(d8_list, dtype=np.int32)
    ei_patch = {}
    for (rs, cs) in bfs_extent:
        if (rs, cs) in d8_ei:
            ei_patch[(rs, cs)] = d8_ei[(rs, cs)]
        else:
            dists   = (d8_arr[:, 0] - rs) ** 2 + (d8_arr[:, 1] - cs) ** 2
            nearest = d8_list[int(np.argmin(dists))]
            ei_patch[(rs, cs)] = d8_ei[nearest]
    return ei_patch
def _run_aoi_model(name, dem_arr, flow_arr, aoi_arr, nrows, ncols,
                   pixel_size, tana, bfs_max_n):
    window     = RECOVERY_WINDOW_M.get(name, 0.0)
    ei_arr     = np.zeros((nrows, ncols), dtype=np.float32)
    candidates = np.transpose(np.nonzero((aoi_arr >= 1) & np.isfinite(dem_arr)))
    arcpy.AddMessage("    {} candidate pixels.".format(len(candidates)))
    for idx, (r0, c0) in enumerate(candidates):
        if idx % 500 == 0:
            arcpy.AddMessage("    {} {}/{}...".format(name, idx, len(candidates)))
        if not np.isfinite(dem_arr[r0, c0]):
            continue
        for (r, c), v in sim_path(r0, c0, dem_arr[r0, c0], tana,
                                   dem_arr, flow_arr, nrows, ncols,
                                   pixel_size, bfs_max_n,
                                   recovery_window=window).items():
            ei_arr[r, c] += v
    return ei_arr
def run_ria(dem_arr, flow_arr, prone_arr, nrows, ncols, pixel_size,
            tana_ria, min_patch_ria, bfs_max_n):
    arcpy.AddMessage("  Running RIA...")
    clean = prepare_ria_aoi(prone_arr, min_patch_ria)
    ei    = _run_aoi_model("RIA", dem_arr, flow_arr, clean, nrows, ncols,
                           pixel_size, tana_ria, bfs_max_n)
    arcpy.AddMessage("    RIA done.")
    return ei
def run_glof(dem_arr, flow_arr, lake_arr, nrows, ncols, pixel_size,
             tana_glof, bfs_max_n):
    arcpy.AddMessage("  Running GLOF...")
    window     = RECOVERY_WINDOW_M.get("GLOF", 0.0)
    struct     = generate_binary_structure(2, 2)
    lbl, _     = label((lake_arr > 0).astype(np.uint8), structure=struct)
    unique_ids = np.unique(lbl[lbl > 0])
    ei_arr     = np.zeros((nrows, ncols), dtype=np.float32)
    arcpy.AddMessage("    {} lakes found.".format(len(unique_ids)))
    for count, lid in enumerate(unique_ids, 1):
        if count % 10 == 0 or count == len(unique_ids):
            arcpy.AddMessage("    GLOF lake {}/{}...".format(count, len(unique_ids)))
        mask      = (lbl == lid)
        lake_size = float(np.sum(mask))
        search = (binary_dilation(mask, structure=struct, iterations=GLOF_DILATION) & ~mask
                  if GLOF_DILATION > 0 else mask)
        valid  = [(dem_arr[r, c], r, c) for r, c in np.argwhere(search)
                  if 0 <= r < nrows and 0 <= c < ncols and np.isfinite(dem_arr[r, c])]
        if not valid:
            continue
        elev0, r0, c0 = min(valid, key=lambda x: x[0])
        weight = math.sqrt(math.log(lake_size + 1.0))
        for (r, c), v in sim_path(r0, c0, elev0, tana_glof,
                                   dem_arr, flow_arr, nrows, ncols,
                                   pixel_size, bfs_max_n, weight=weight,
                                   recovery_window=window).items():
            ei_arr[r, c] += v
    arcpy.AddMessage("    GLOF done.")
    return ei_arr
def run_ls(dem_arr, flow_arr, ls_arr, ria_arr, nrows, ncols, pixel_size,
           tana_ls, min_patch_ls, bfs_max_n):
    arcpy.AddMessage("  Running LS...")
    clean = prepare_ls_aoi(ls_arr, ria_arr, min_patch_ls)
    ei    = _run_aoi_model("LS", dem_arr, flow_arr, clean, nrows, ncols,
                           pixel_size, tana_ls, bfs_max_n)
    arcpy.AddMessage("    LS done.")
    return ei
FIELD_MAP = {"RIA": "EI_RIA", "GLOF": "EI_GLOF", "LS": "EI_LS", "Multi": "EI_Multi"}
CLASS_FIELD_MAP = {"RIA": "CL_RIA", "GLOF": "CL_GLOF", "LS": "CL_LS", "Multi": "CL_Multi"}
def classify_ei(ei_value):
    if ei_value <= 0:
        return 0
    if ei_value <= 0.2:
        return 1
    if ei_value <= 0.4:
        return 2
    if ei_value <= 0.6:
        return 3
    if ei_value <= 0.8:
        return 4
    return 5
def make_transform_inv(x_min, y_top, pixel_size):
    def _inv(x, y):
        return (x - x_min) / pixel_size, (y_top - y) / pixel_size
    return _inv
def _sample_hi_at_point(x, y, hi_arr, transform_inv, nrows, ncols):
    col_f, row_f = transform_inv(x, y)
    row_i, col_i = int(math.floor(row_f)), int(math.floor(col_f))
    if 0 <= row_i < nrows and 0 <= col_i < ncols:
        v = float(hi_arr[row_i, col_i])
        return v if v > 0 else -1.0
    return -1.0
def _densify_polyline(geometry, densify_dist):
    densified = geometry.densify("DISTANCE", densify_dist)
    pts = []
    for part in densified:
        for pt in part:
            if pt:
                pts.append((pt.X, pt.Y))
    return pts
def _ei_for_polyline(geometry, hi_arr, transform_inv, nrows, ncols, densify_dist):
    pts = _densify_polyline(geometry, densify_dist)
    if not pts:
        return -1.0, -1.0
    vals = [_sample_hi_at_point(x, y, hi_arr, transform_inv, nrows, ncols)
            for x, y in pts]
    positive = [v for v in vals if v > 0]
    if not positive:
        return -1.0, -1.0
    # spacing is uniform, so the simple mean is the length-weighted mean
    return float(np.mean(positive)), float(np.max(positive))
def _ei_for_polygon(geometry, hi_arr, transform_inv, nrows, ncols):
    c = geometry.centroid
    return _sample_hi_at_point(c.X, c.Y, hi_arr, transform_inv, nrows, ncols)
def run_exposure(infra_paths, hi_dict, ll, pixel_size, nrows, ncols,
                 out_folder, sr):
    if not infra_paths:
        return
    arcpy.AddMessage("PART 5: Exposure assessment...")
    transform_inv = make_transform_inv(ll.X, ll.Y + nrows * pixel_size, pixel_size)
    densify_dist  = pixel_size
    active_tags   = list(hi_dict.keys())
    for infra_path in infra_paths:
        layer_name = os.path.splitext(os.path.basename(infra_path))[0]
        arcpy.AddMessage("  Processing: {}".format(layer_name))
        desc      = arcpy.Describe(infra_path)
        geom_type = desc.shapeType
        src_sr    = desc.spatialReference
        if src_sr.factoryCode != sr.factoryCode:
            arcpy.AddMessage("    Reprojecting {} to HI raster CRS...".format(layer_name))
            reproj_path = os.path.join(out_folder, layer_name + "_reproj.shp")
            arcpy.Project_management(infra_path, reproj_path, sr)
            work_path = reproj_path
        else:
            work_path = infra_path
        out_path = os.path.join(out_folder, layer_name + "_exposure.shp")
        arcpy.CopyFeatures_management(work_path, out_path)
        for tag in active_tags:
            field_name = FIELD_MAP[tag]
            arcpy.AddField_management(out_path, field_name, "FLOAT")
            if geom_type == "Polyline":
                arcpy.AddField_management(out_path, field_name + "_mx", "FLOAT")
            arcpy.AddField_management(out_path, CLASS_FIELD_MAP[tag], "SHORT")
        arcpy.AddMessage("    Geometry type: {} | EI fields: {} | Class fields: {}".format(
            geom_type,
            ", ".join(FIELD_MAP[t] for t in active_tags),
            ", ".join(CLASS_FIELD_MAP[t] for t in active_tags)))
        update_fields = ["SHAPE@"]
        for tag in active_tags:
            update_fields.append(FIELD_MAP[tag])
            if geom_type == "Polyline":
                update_fields.append(FIELD_MAP[tag] + "_mx")
            update_fields.append(CLASS_FIELD_MAP[tag])
        n_features = int(arcpy.GetCount_management(out_path)[0])
        arcpy.AddMessage("    {} features to process.".format(n_features))
        with arcpy.da.UpdateCursor(out_path, update_fields) as cursor:
            for feat_idx, row in enumerate(cursor):
                if feat_idx % 200 == 0:
                    arcpy.AddMessage("    {}/{} features...".format(feat_idx, n_features))
                geom = row[0]
                if geom is None:
                    for i in range(1, len(update_fields)):
                        row[i] = 0 if update_fields[i] in CLASS_FIELD_MAP.values() else -1.0
                    cursor.updateRow(row)
                    continue
                field_idx = 1
                for tag in active_tags:
                    hi_arr = hi_dict[tag]
                    if geom_type == "Point":
                        pt  = geom.firstPoint
                        val = _sample_hi_at_point(pt.X, pt.Y, hi_arr,
                                                  transform_inv, nrows, ncols)
                        row[field_idx]     = val
                        row[field_idx + 1] = classify_ei(val)
                        field_idx += 2
                    elif geom_type == "Polyline":
                        mean_ei, max_ei = _ei_for_polyline(geom, hi_arr, transform_inv,
                                                           nrows, ncols, densify_dist)
                        row[field_idx]     = mean_ei
                        row[field_idx + 1] = max_ei
                        row[field_idx + 2] = classify_ei(mean_ei)
                        field_idx += 3
                    elif geom_type == "Polygon":
                        val = _ei_for_polygon(geom, hi_arr, transform_inv, nrows, ncols)
                        row[field_idx]     = val
                        row[field_idx + 1] = classify_ei(val)
                        field_idx += 2
                    else:
                        row[field_idx]     = -1.0
                        row[field_idx + 1] = 0
                        field_idx += 2
                cursor.updateRow(row)
        arcpy.AddMessage("    Saved: {}".format(os.path.basename(out_path)))
    arcpy.AddMessage("  Exposure assessment complete.")
def main():
    input_dem  = arcpy.GetParameterAsText(0)
    input_ria  = arcpy.GetParameterAsText(1).strip() or None
    input_glof = arcpy.GetParameterAsText(2).strip() or None
    input_ls   = arcpy.GetParameterAsText(3).strip() or None
    out_folder = arcpy.GetParameterAsText(4)
    active_names = [h for h, v in
                    [("RIA", input_ria), ("GLOF", input_glof), ("LS", input_ls)] if v]
    if not active_names:
        arcpy.AddError("At least one AOI input (RIA, GLOF, or LS) is required.")
        sys.exit(1)
    arcpy.AddMessage("Active hazards: {}".format(", ".join(active_names)))
    def parse_float(idx, default):
        try:
            v = arcpy.GetParameterAsText(idx).strip()
            return float(v) if v else default
        except (ValueError, AttributeError):
            return default
    def parse_int(idx, default, min_val=1):
        try:
            v = arcpy.GetParameterAsText(idx).strip()
            return max(min_val, int(v)) if v else default
        except (ValueError, AttributeError):
            return default
    tana_ria      = parse_float(5, DEFAULT_TANA_RIA)
    tana_glof     = parse_float(6, DEFAULT_TANA_GLOF)
    tana_ls       = parse_float(7, DEFAULT_TANA_LS)
    bfs_max_n     = parse_int(8,  DEFAULT_BFS_MAX_N)
    min_patch_ria = parse_int(9,  DEFAULT_MIN_PATCH_RIA, min_val=1)
    min_patch_ls  = parse_int(10, DEFAULT_MIN_PATCH_LS,  min_val=1)
    infra_raw   = arcpy.GetParameterAsText(11).strip()
    infra_paths = [p.strip().strip("'\"") for p in infra_raw.split(";")
                   if p.strip().strip("'\"")] if infra_raw else []
    if infra_paths:
        arcpy.AddMessage("Infrastructure layers: {}".format(
            ", ".join(os.path.basename(p) for p in infra_paths)))
    else:
        arcpy.AddMessage("No infrastructure layers provided - skipping exposure assessment.")
    arcpy.env.overwriteOutput = True
    temp_folder = os.path.join(out_folder, "MH_temp")
    os.makedirs(temp_folder, exist_ok=True)
    os.makedirs(out_folder,  exist_ok=True)
    def tmp(n):  return os.path.join(temp_folder, n)
    def outp(n): return os.path.join(out_folder, n)
    try:
        dem_reproj, dem_reason = needs_projection(input_dem)
        if dem_reproj:
            arcpy.AddMessage("Reprojecting DEM ({})".format(dem_reason))
            clon, clat = get_raster_center(input_dem)
            utm_sr, _  = build_utm_sr(clon, clat)
            dem_work   = tmp("dem_utm.tif")
            arcpy.ProjectRaster_management(input_dem, dem_work, utm_sr, "BILINEAR")
        else:
            dem_work = input_dem
        ria_work  = align_raster(input_ria,  dem_work, tmp("ria_reproj.tif"))  if input_ria  else None
        glof_work = align_raster(input_glof, dem_work, tmp("glof_reproj.tif")) if input_glof else None
        ls_work   = align_raster(input_ls,   dem_work, tmp("ls_reproj.tif"))   if input_ls   else None
        arcpy.env.snapRaster = dem_work
        arcpy.env.cellSize   = dem_work
        arcpy.env.extent     = dem_work
        arcpy.env.mask       = dem_work
        arcpy.AddMessage("PART 1: Fill -> FocalStatistics -> Fill -> FlowDirection...")
        arcpy.AddMessage("  Recovery window per hazard: {}".format(
            {k: ("unbounded" if v is None else "{:.0f} m".format(v))
             for k, v in RECOVERY_WINDOW_M.items()}))
        dem_fill1  = Fill(Raster(dem_work))
        dem_smooth = dem_fill1
        for i in range(FOCAL_ITERATIONS):
            dem_smooth = FocalStatistics(dem_smooth, NbrRectangle(3, 3, "CELL"), "MEAN")
            arcpy.AddMessage("  FocalStatistics pass {}/{} done.".format(i + 1, FOCAL_ITERATIONS))
        dem_cr = Fill(dem_smooth)
        dem_cr.save(tmp("dem_cr.tif"))
        fldir = FlowDirection(dem_cr)
        fldir.save(tmp("fldir.tif"))
        pixel_size = float(arcpy.GetRasterProperties_management(
            tmp("dem_cr.tif"), "CELLSIZEX").getOutput(0))
        arcpy.AddMessage("  Pixel size: {} m".format(pixel_size))
        arcpy.AddMessage("PART 2: Reading arrays...")
        dem_cr_r = arcpy.Raster(tmp("dem_cr.tif"))
        nrows    = dem_cr_r.height
        ncols    = dem_cr_r.width
        ll       = arcpy.Point(dem_cr_r.extent.XMin, dem_cr_r.extent.YMin)
        sr       = dem_cr_r.spatialReference
        dem_arr  = arcpy.RasterToNumPyArray(dem_cr_r, ll, ncols, nrows,
                                            nodata_to_value=-9999).astype(np.float32)
        dem_arr[dem_arr == -9999] = np.nan
        flow_arr = arcpy.RasterToNumPyArray(tmp("fldir.tif"), ll, ncols, nrows,
                                            nodata_to_value=0).astype(np.int32)
        def read_prone(work_path, bin_name):
            Con(IsNull(Raster(work_path)), 0, 1).save(tmp(bin_name))
            return arcpy.RasterToNumPyArray(tmp(bin_name), ll, ncols, nrows,
                                            nodata_to_value=0).astype(np.int32)
        ria_arr  = read_prone(ria_work,  "ria_bin.tif")  if ria_work  else None
        glof_arr = read_prone(glof_work, "glof_bin.tif") if glof_work else None
        ls_arr   = read_prone(ls_work,   "ls_bin.tif")   if ls_work   else None
        arcpy.AddMessage("  Grid: {}x{}".format(nrows, ncols))
        arcpy.AddMessage("PART 3: Running hazard models...")
        ei_ria  = run_ria (dem_arr, flow_arr, ria_arr,  nrows, ncols, pixel_size,
                           tana_ria,  min_patch_ria, bfs_max_n) if ria_arr  is not None else None
        ei_glof = run_glof(dem_arr, flow_arr, glof_arr, nrows, ncols, pixel_size,
                           tana_glof, bfs_max_n)                if glof_arr is not None else None
        ei_ls   = run_ls  (dem_arr, flow_arr, ls_arr,   ria_arr, nrows, ncols, pixel_size,
                           tana_ls,   min_patch_ls,  bfs_max_n) if ls_arr   is not None else None
        arcpy.AddMessage("PART 4: Saving HI raster outputs...")
        hi_norm = {}
        ei_components = []
        for ei, tag, fname in [(ei_ria,  "RIA",  "RIA_HI.tif"),
                               (ei_glof, "GLOF", "GLOF_HI.tif"),
                               (ei_ls,   "LS",   "LS_HI.tif")]:
            if ei is None:
                continue
            raw_name = fname.replace("_HI.tif", "_HI_raw.tif")
            save_raster(ei.astype(np.float32), ll, pixel_size, sr, outp(raw_name))
            arcpy.AddMessage("  {} saved.".format(raw_name))
            ei_n = normalize(ei)
            save_raster(ei_n, ll, pixel_size, sr, outp(fname))
            arcpy.AddMessage("  {} saved.".format(fname))
            hi_norm[tag] = ei_n
            ei_components.append(ei_n)
        if len(ei_components) >= 2:
            multi_n = normalize(sum(ei_components))
            save_raster(multi_n, ll, pixel_size, sr, outp("Multi_HI.tif"))
            arcpy.AddMessage("  Multi_HI.tif saved.")
            hi_norm["Multi"] = multi_n
        else:
            arcpy.AddMessage("  Only one hazard active - Multi_HI not generated.")
        if infra_paths and hi_norm:
            run_exposure(infra_paths, hi_norm, ll, pixel_size, nrows, ncols,
                         out_folder, sr)
        elif infra_paths and not hi_norm:
            arcpy.AddWarning("Infrastructure layers provided but no HI rasters were "
                             "produced. Check that at least one AOI input is valid.")
        arcpy.AddMessage("Done. Outputs: {}".format(out_folder))
        arcpy.AddMessage(
            "  Active: {} | tan alpha RIA={} GLOF={} LS={} | BFS max N={} | "
            "Min patch RIA={} px LS={} px".format(
                ", ".join(active_names), tana_ria, tana_glof, tana_ls,
                bfs_max_n, min_patch_ria, min_patch_ls))
        if infra_paths:
            arcpy.AddMessage("  Exposure outputs: {} layer(s)".format(len(infra_paths)))
            arcpy.AddMessage("  EI fields: {} | Class fields: {}".format(
                ", ".join(FIELD_MAP[t] for t in hi_norm),
                ", ".join(CLASS_FIELD_MAP[t] for t in hi_norm)))
    except arcpy.ExecuteError:
        arcpy.AddError(arcpy.GetMessages(2))
        sys.exit(1)
    except Exception:
        arcpy.AddError(traceback.format_exc())
        sys.exit(1)
    finally:
        arcpy.CheckInExtension("Spatial")
if __name__ == "__main__":
    main()
