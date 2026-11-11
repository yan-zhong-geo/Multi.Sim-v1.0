<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/multisim-logo-dark.png">
    <img src="assets/Multi.Sim Logo.png" alt="Multi.Sim logo" width="600">
  </picture>
</p>

# Multi.Sim v1.0

**A GIS toolbox for regional runout screening and exposure assessment**

**Developer:** [Yan Zhong](https://sites.google.com/view/yanzhong-geo)  
**Affiliation:** [Climate Change Impacts and Risks in the Anthropocene (C-CIA), University of Geneva](https://c-cia.ch/)  
**Contact:** yan.zhong@unige.ch; yan.zhong.geo@gmail.com  
**Platform:** ArcGIS Pro Python toolbox (`.atbx`)  
**Toolbox:** `Multi.Sim_v1.0.atbx`

---

## 1. What Multi.Sim does

Multi.Sim v1.0 is a regional-scale GIS tool for screening potential runout corridors associated with three types of mass-movement hazards:

- **Rock-ice avalanches (RIA)**
- **Glacial lake outburst floods (GLOF)**
- **Long-runout landslides (LS)**, including flow-type events such as debris flows

It uses a DEM, mapped source areas or lake masks, and a simplified terrain-routing approach to generate a **Hazard Index (HI)** raster for each selected hazard. When two or more hazards are simulated, it also creates a composite **Multi-Hazard HI** raster. The toolbox can then sample these rasters against infrastructure or other exposed features to calculate an **Exposure Index (EI)**.

Multi.Sim is intended for **first-pass regional screening**: identifying potentially affected valleys, corridors, assets, and locations that may need field checking or more detailed modelling. It is not an event-scale physical model and should not be used on its own for engineering design, evacuation timing, or definitive risk estimates.

> **HI and EI are relative indices.** They do not represent flow depth, velocity, discharge, impact pressure, probability, or risk. A high value indicates a relatively high position within the modelled index surface, not a quantified physical impact.

## 2. Toolbox contents and which tool to use

The current toolbox contains **two tools**:

| Tool | Use it when… | Main outputs |
|---|---|---|
| **1. Multi.Sim v1.0** | You need to generate hazard runout HI rasters from a DEM and source data. | Raw and normalized HI rasters; optional feature-level exposure outputs. |
| **2. Exposure Assessment** | You already have HI rasters and want to assess new exposure layers or summarize exposure by reporting zones. | Feature-level EI outputs; optional aggregation summary. |

A typical project uses the first tool once to generate the HI rasters, then uses the second tool whenever you want to assess another set of buildings, roads, bridges, or other assets without rerunning the hazard simulations.

## 3. Requirements

- ArcGIS Pro
- Spatial Analyst extension
- Python package: `scipy`
- A DEM covering the study area
- At least one hazard source raster for **Multi.Sim v1.0**
- Optional point, polyline, or polygon exposure layers

The toolbox processes terrain distances in metres. It reprojects geographic or non-metre projected DEMs to an appropriate UTM coordinate system. **If a dataset has an unknown or incorrectly defined coordinate system, define or repair its coordinate system before running the tool.** Reprojection cannot correct an incorrectly assigned CRS.

## 4. Before you run: prepare your data

Good input preparation is essential. Check the following before opening the toolbox.

### 4.1 DEM

1. Use a DEM that covers the full source-to-runout area, including the downstream valley or fan where runout may extend.
2. Check that the DEM has the correct coordinate system, vertical units, and elevation values.
3. Inspect for NoData gaps, striping, artificial pits, abrupt tile boundaries, and other artefacts that may disrupt flow routing.
4. Use a consistent DEM for all hazards in the same study so their outputs share the same terrain basis and grid.

The tool conditions the DEM by filling sinks, applying one 3 × 3 focal-mean smoothing pass, filling sinks again, and deriving D8 flow directions. This helps reduce small terrain artefacts, but it cannot fix major errors in the input DEM.

### 4.2 RIA and landslide source rasters

For RIA and LS, **every valid (non-NoData) raster cell is treated as a source cell**. The cell value itself is not used as a probability or intensity weight.

- Prepare a source mask containing valid cells only where you want the model to start.
- Set cells outside the source area to **NoData**, not zero, unless zero is intentionally part of the source mask.
- If your starting dataset is a continuous susceptibility map, choose and apply a threshold before running Multi.Sim, then convert cells outside the selected source area to NoData.
- Review the resulting source mask on top of the DEM. An overly broad source mask can create a very large number of seeds and extensive runout corridors.

Small contiguous source patches are filtered using the minimum-patch parameters (defaults: 22 pixels for RIA and 4 pixels for LS). These are pixel counts, not map areas.

### 4.3 GLOF lake raster

For GLOF, valid raster cells represent lake areas. The tool identifies individual lakes using **8-connectivity**.

- Ensure each lake is represented by a contiguous group of valid cells.
- Under 8-connectivity, cells touching diagonally belong to the same lake. If two lakes should be separate, check that their masks do not touch diagonally.
- The tool searches around each lake boundary for a low-elevation outlet. The default boundary search-ring dilation is 5 pixels.
- Check the automatically identified outlets against the DEM and mapped drainage routes, especially where lake boundaries are noisy or the DEM has artefacts.

### 4.4 Exposure layers

Exposure layers can be points, polylines, or polygons. Examples include building footprints, road or trail networks, bridges, power lines, and other infrastructure.

- Make sure geometries are valid and the layer has a correctly defined coordinate system.
- Use a layer that represents the assets you actually want to assess; remove duplicates and irrelevant features first.
- The tool samples points at their locations, polygons at their centroids, and lines along their length. These choices matter when features are large or spatially complex (see [Section 7](#7-how-exposure-is-calculated)).

## 5. Workflow 1 — Multi.Sim v1.0: generate hazard maps

Use **1. Multi.Sim v1.0** when you have a DEM and at least one hazard source raster and want to generate new HI maps.

### 5.1 Parameters

| Parameter | Required? | Description / default |
|---|---:|---|
| `DEM` | Yes | DEM used for terrain conditioning and flow routing. |
| `RIA_AOI` | No | Raster mask of RIA source cells. |
| `GLOF_AOI` | No | Raster mask of glacial lake areas. |
| `LS_AOI` | No | Raster mask of landslide source cells. |
| `Output Folder` | Yes | Folder for output rasters, feature layers, and intermediate files. |
| `Tana RIA` | No | RIA reach-angle tangent; default `0.10` (about 6°). |
| `Tana GLOF` | No | GLOF reach-angle tangent; default `0.05` (about 3°). |
| `Tana LS` | No | LS reach-angle tangent; default `0.19` (about 11°). |
| `BFS max neighbours` | No | Maximum number of steepest downslope neighbours used for lateral spread; default `2`. |
| `Min patch pixels RIA` | No | Minimum contiguous RIA source-patch size; default `22` pixels. |
| `Min patch pixels LS` | No | Minimum contiguous LS source-patch size; default `4` pixels. |
| `Exposure` | No | Optional vector feature layer for exposure assessment in the same run. |

At least one of `RIA_AOI`, `GLOF_AOI`, or `LS_AOI` must be supplied. Only the selected hazard models are run.

### 5.2 Running the tool

1. Add `Multi.Sim_v1.0.atbx` to your ArcGIS Pro project.
2. Open **1. Multi.Sim v1.0** in the toolbox.
3. Select the DEM and one or more prepared hazard source rasters.
4. Choose an output folder with enough free disk space.
5. Start with the default parameters unless you have a reason to change them. See [Parameter guidance](#8-parameter-guidance) before adjusting them.
6. Optionally add an exposure feature layer if you want EI fields calculated in the same run.
7. Run the tool. Review the geoprocessing messages and confirm that the expected outputs have been written.
8. Add the normalized HI rasters to the map and inspect them against the DEM, source areas, rivers, and known terrain constraints. Do not rely on the raster statistics alone.

### 5.3 What the model does

The model uses a shared three-part procedure for the active hazards:

**1. Trace a downslope D8 path.** Each RIA or LS source cell is processed as a seed. For GLOF, the seed is the lowest-elevation outlet cell identified around each lake. The model follows the DEM-derived D8 flow direction and evaluates the cumulative reach angle against the selected `Tana` threshold.

**2. Expand the corridor laterally.** A breadth-first search (BFS) spreads from the D8 path to up to the configured number of steepest downslope neighbours. Candidate cells must be at the same or a lower elevation. This creates a simplified runout corridor around the main path.

**3. Assign and accumulate HI.** Cells in the corridor receive a value based on their nearest D8-path position. The value decreases quadratically with distance along the path, and contributions from multiple seeds are accumulated. The accumulated surface is then normalized for mapping.

The D8 path does not always stop at the first local reach-angle threshold violation. A recovery window allows the path to continue through low-gradient sections before termination. The current windows are 3 km for RIA and LS and unbounded for GLOF.

### 5.4 Hazard-specific behaviour

| Hazard | Seed definition | Additional behaviour |
|---|---|---|
| RIA | Every valid RIA source cell | Small source patches are removed; contributions from source cells are accumulated. |
| GLOF | Lowest outlet cell around each 8-connected lake | Each lake is weighted by the square root of the natural logarithm of its cell count plus one; lake contributions are accumulated. |
| LS | Every valid LS source cell | Small source patches are removed; cells overlapping RIA source areas are excluded. The model is intended for long-runout, flow-type landslides, not every landslide mechanism. |

### 5.5 Outputs

Only outputs for active hazards are created.

| Output | Meaning |
|---|---|
| `RIA_HI_raw.tif` | Accumulated RIA index before normalization. |
| `RIA_HI.tif` | Normalized RIA HI (0–1). |
| `GLOF_HI_raw.tif` | Accumulated GLOF index before normalization. |
| `GLOF_HI.tif` | Normalized GLOF HI (0–1). |
| `LS_HI_raw.tif` | Accumulated LS index before normalization. |
| `LS_HI.tif` | Normalized LS HI (0–1). |
| `Multi_HI.tif` | Composite Multi-Hazard HI, created when at least two hazards are active. |
| `MH_temp/` | Intermediate processing files. |

The normalized HI is produced by clipping positive accumulated values at the 99th percentile and rescaling the result to 0–1. The composite `Multi_HI.tif` is made by adding the normalized individual hazard rasters and normalizing the sum again.

**Important:** `Multi_HI.tif` is an additive spatial overlay. It does not simulate physical triggering or interactions between hazards (for example, an avalanche entering a lake and causing an outburst flood).

### 5.6 Optional: exposure in the same run

If you provide an `Exposure` layer, the tool creates an output feature layer named `<LayerName>_exposure.shp` and adds EI and class fields for each active hazard (and Multi, if available):

- `EI_RIA`, `EI_GLOF`, `EI_LS`, `EI_Multi`: element-level EI values.
- `CL_RIA`, `CL_GLOF`, `CL_LS`, `CL_Multi`: exposure classes (0–5).
- `EI_*_mx`: maximum EI along polyline features, where applicable.

For details on the sampling rules and classes, see [How exposure is calculated](#7-how-exposure-is-calculated).

## 6. Workflow 2 — Exposure Assessment: assess existing hazard maps

Use **2. Exposure Assessment** if you already have one or more HI rasters and want to assess additional assets or summarize exposure by administrative units, watersheds, or other polygons. This tool does **not** rerun the hazard simulations.

### 6.1 Inputs

| Parameter | Description |
|---|---|
| `RIA HI Raster` | Optional normalized RIA HI raster. |
| `GLOF HI Raster` | Optional normalized GLOF HI raster. |
| `LS HI Raster` | Optional normalized LS HI raster. |
| `Multi HI Raster` | Optional composite Multi-Hazard HI raster. |
| `Exposure Layers` | Optional one or more point, polyline, or polygon feature layers. |
| `Aggregation Layer` | Optional polygon layer defining reporting units. |
| `Output Folder` | Required output folder. |

At least one HI raster is required. Supply exposure layers to create feature-level EI outputs. Supply an aggregation layer if you also want a summary by reporting unit.

### 6.2 Running the tool

1. Open **2. Exposure Assessment**.
2. Select the normalized HI raster(s) you want to use. Do not select raw accumulated HI rasters unless you specifically intend to assess those non-normalized values.
3. Add one or more exposure layers.
4. If required, add a polygon aggregation layer (for example, municipalities or sub-basins).
5. Choose an output folder and run the tool.
6. Inspect the output feature layers and, if requested, the aggregation summary. Check that the input HI rasters and exposure data overlap and use the intended study area.

You can rerun this workflow with different exposure datasets while reusing the same HI rasters. This is useful when maintaining separate inventories for buildings, roads, bridges, and other asset types.

### 6.3 Outputs

For each processed exposure layer, the tool creates `<LayerName>_exposure.shp` with EI and class fields. For shapefile compatibility, field names are shortened:

- `EIRIA`, `EIGLOF`, `EILS`, `EIMUL`: EI values for the corresponding HI raster.
- `CL_RIA`, `CL_GLOF`, `CL_LS`, `CL_Multi`: exposure classes (0–5).
- `EIRIAMX`, `EIGLOFMX`, `EILSMX`, `EIMULMX`: maximum EI along polylines, where applicable.

If an aggregation layer is supplied, the tool also creates:

- `aggregation_summary.shp`: aggregation polygons with summed positive EI values for the processed exposure layers and hazards;
- `field_legend.csv`: descriptions of the summary fields;
- `EXP_temp/`: intermediate aggregation files.

The aggregation summary is a sum of positive element-level EI values within each reporting polygon. It is useful for comparing the relative concentration of mapped exposure in the selected units, but it is not a population count, expected loss, or risk estimate. The result depends on the number and segmentation of features in the input inventory.

## 7. How exposure is calculated

EI is obtained by sampling the selected HI raster(s) against each vector feature:

- **Point:** sample the raster cell at the point location.
- **Polygon:** sample the raster at the polygon centroid. This is not an average over the full polygon.
- **Polyline:** densify the line at approximately one raster-cell spacing, then calculate the mean and maximum sampled HI. The exposure class is based on the mean EI.

A feature with no positive sampled value receives the no-exposure sentinel value `-1` in the EI field and class `0`. Check this when filtering or symbolizing the output: `-1` is a no-positive-sample marker, not a negative exposure value.

### Exposure classes

| Class | Label | EI range |
|---:|---|---|
| 0 | Negligible / no positive sample | EI = 0 or no positive sampled value |
| 1 | Very low | 0 < EI ≤ 0.2 |
| 2 | Low | 0.2 < EI ≤ 0.4 |
| 3 | Moderate | 0.4 < EI ≤ 0.6 |
| 4 | High | 0.6 < EI ≤ 0.8 |
| 5 | Very high | EI > 0.8 |

These are relative ranking classes applied to normalized HI values. They are not physical impact thresholds or probabilities. Polygon-centroid sampling can miss hazard overlap elsewhere in a large polygon; for line features, the mean and maximum provide complementary information, so consider both when reviewing critical infrastructure.

## 8. Parameter guidance

### 8.1 `Tana` (Fahrböschung threshold)

`Tana` is the tangent of the reach angle used to assess whether a path can continue downslope. Lower values generally allow longer runout; higher values generally shorten it.

| Hazard | Default `Tana` | Approximate angle |
|---|---:|---:|
| RIA | 0.10 | 6° |
| GLOF | 0.05 | 3° |
| LS | 0.19 | 11° |

Treat these as starting values, not universal constants. If you have a regional inventory of observed runout distances, test and calibrate the parameter against those observations. Avoid tuning a value solely to make one desired event footprint match.

### 8.2 `BFS max neighbours`

The default is `2`. Increasing this value generally expands the lateral corridor by allowing more downslope neighbours to be explored. It does not increase the D8 main-path runout distance. Because the search only accepts equal- or lower-elevation cells, it cannot represent volume-driven overtopping of higher valley margins or super-elevation.

### 8.3 Minimum source-patch size

Defaults are 22 pixels for RIA and 4 pixels for LS. These thresholds were developed for a 30 m grid and are based on pixel counts. If you change raster resolution, reconsider the thresholds in terms of the physical source-area size you intend to retain. GLOF lakes are not filtered by these minimum-patch parameters.

### 8.4 Raster resolution and alignment

The default pixel-based settings were developed around a 30 m DEM. Finer or coarser grids change the physical meaning of a pixel-based BFS width and minimum patch size. Keep the DEM and source masks aligned, and review the resampled source masks after processing. Do not assume that changing the cell size leaves the model behaviour unchanged.

## 9. Quality checks and common pitfalls

Before interpreting the maps, check the following:

- **CRS:** all inputs have correctly defined coordinate systems. An unknown or incorrectly assigned CRS should be fixed before processing.
- **Extent:** the DEM covers the entire possible runout path, not just the source area.
- **NoData:** RIA and LS source cells are valid pixels; outside-source cells are NoData. Do not use a background value that the tool will interpret as a source.
- **Lake connectivity:** diagonally touching lake masks are treated as one lake under 8-connectivity.
- **Outlet locations:** inspect GLOF outlet cells, particularly for complex lake boundaries or low-quality DEMs.
- **Exposure geometry:** polygons are sampled at centroids, and polylines at roughly one-cell intervals. Confirm that this is appropriate for the feature type and raster resolution.
- **Output interpretation:** normalized HI and EI are relative indices. Do not label them as inundation depth, velocity, probability, or risk.
- **Cross-study comparison:** each normalized raster is scaled using its own positive-value distribution. Values from different regions or runs are not automatically directly comparable.

## 10. Model scope and limitations

- The model is designed for regional screening, not detailed event reconstruction.
- RIA, GLOF, and LS are simulated independently. Cascading triggers and physical interactions are not explicitly represented.
- The lateral spread is a simplified deterministic BFS and is limited to equal- or lower-elevation neighbours.
- The tool does not calculate flow depth, velocity, discharge, travel time, impact pressure, or warning time.
- Results depend on DEM quality, source-area delineation, lake inventories, raster resolution, and exposure-data completeness and positional accuracy.
- The LS parameterization is intended for long-runout, flow-type landslides and should not be assumed to represent all landslide mechanisms.

Use the outputs to prioritize field checks, data improvement, and more detailed site-specific modelling.

## 11. Publication

Please cite the accompanying paper when using Multi.Sim:

> Zhong, Y., Allen, S., Steiner, J., Kincey, M., Rosser, N., Munch, J., Adhikari, B. R., Pandey, V. P., & Stoffel, M. (2026). *Multi.Sim v1.0: A Multi-Hazard Simulation Framework for Runout Modeling and Exposure Assessment in High Mountain Regions.*

## 12. License

See [`LICENSE`](LICENSE) for the terms of use.

