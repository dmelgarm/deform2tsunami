# deform2tsunami

Export triangular surface deformation and tsunami initial conditions to WGS84 longitude–latitude GeoTIFFs. The command and Python function use the same sampler. Inputs use the PyLith-style surface layout described below.

## Install

From a downloaded or cloned copy, in a Python environment of your choice:

```sh
python -m pip install .
deform2tsunami --help
```

The four runtime dependencies are NumPy, h5py, Matplotlib (triangle lookup), and Rasterio (coordinates and GeoTIFF). Rasterio wheels include GDAL for common platforms. Jupyter is needed only to run the [example notebook](extract_rasters.ipynb).

## Export a raster

```sh
deform2tsunami raster surface.h5 --field vertical --out uplift.tif
deform2tsunami raster surface.h5 --field bed-change --out bed_change.tif
deform2tsunami raster surface.h5 --field eta0 --out initial_surface.tif
deform2tsunami kajiura surface.h5 --bounds -125.8 44.5 -124.1 45.6 --arcsec 30 --out recomputed.tif
```

| Command / field | Value in metres, positive up |
| --- | --- |
| `raster --field vertical` | Vertical ground displacement on land and offshore. |
| `raster --field bed-change` | First-order change in bed or land elevation at a fixed location, including horizontal displacement across the original triangle slope. |
| `raster --field eta0` | Archived initial water-surface elevation interpolated from native vertices. Requires `/tsunami/eta0`. |
| `kajiura` | Initial water-surface elevation recomputed from displacement on a projected 500 m grid, then sampled directly at the requested pixels. Works with a file containing only the original deformation fields. |

Bounds are optional and ordered west, south, east, north in degrees. The default covers the projected mesh bounding rectangle. The default pixel spacing is 30 arcseconds. Raster edges snap outward to multiples of that spacing so identical bounds and spacing give matching rasters across models. Each float32 pixel is a point sample at its center. Rows run north to south. Missing values are NaN. Use `--overwrite` to replace an existing GeoTIFF. Select a regional bounding box for a quick first export.

The same operation in Python:

```python
from deform2tsunami import export_raster

export_raster("surface.h5", "oregon.tif", "bed-change",
              bounds=(-125.8, 44.5, -124.1, 45.6), arcsec=30)
```

Pass `"eta0-recomputed"` as the field for the Python equivalent of the `kajiura` command.

## Input and interpolation

| HDF5 path | Shape | Meaning |
| --- | --- | --- |
| `/geometry/vertices` | (N, 3) | Projected x, projected y, original elevation in metres. |
| `/viz/topology/cells` | (M, 3) | Zero-based indices of the original surface triangles. |
| `/vertex_fields/displacement` | (1, N, 3) | One static displacement per vertex in projected x, projected y and up, in metres. |
| `/tsunami/eta0` | (N,) | Optional archived sea-surface elevations in the same vertex order. NaN denotes dry or unsupported points. |

The horizontal CRS comes from the root `crs_wkt` attribute. For a file without this attribute, supply its known projected metre CRS with `--source-crs` or the Python keyword `source_crs`. This tool expects a triangular, single-valued terrain surface in projected Cartesian metres and does not read a volume mesh, infer a datum from coordinates, or select a frame from a time series.

Each geographic pixel center is projected into the mesh and located on an original triangle. The triangle's original vertex indices preserve displacement jumps at a fault that reaches the surface. A point exactly on a shared edge uses the face selected by the triangle lookup. With bed elevation **b** positive up, bed change is **uz − ux ∂b/∂x − uy ∂b/∂y**. Ground fields cover both wet and dry terrain. Archived water elevations are sampled only on fully submerged triangles with three finite archived values. Mixed wet/dry triangles and unsupported points return NoData. Reducing the output pixel size samples this same piecewise-linear field more densely.

## Kajiura calculation

For rapid bed motion, the constant-depth amplitude transfer is **1/cosh(|k|h)**, with angular wavenumber **k** and water depth **h**. The variable-depth approximation weights each submerged source point between adjacent 250 m reference depths before applying the FFT for each band. Filtering uses a 500 m working grid anchored at the projected origin and a 50 km source margin around each calculation tile. This is a local-depth approximation for tsunami generation. Coastline barriers and propagation require a separate water-domain model.

Dry and out-of-mesh points supply zero bottom forcing. The continuous filtered response is bilinearly sampled at pixel centers whose original terrain elevation is below zero. Around source boundaries the result assumes zero forcing beyond the available mesh. This shoreline convention differs from the conservative fully wet triangle rule for exporting archived values. The working grid and margins are independent of the requested lon/lat pixel spacing. The original permanent bed change and the initial water-surface elevation serve different roles in a tsunami model.

For the homogeneous full-margin T100 example, the notebook covers Cascadia from 130.5°W to 121°W and 39.6°N to 51°N. At 30 arcseconds the rasters have 1,140 columns and 1,368 rows. Approximate maxima are 12.83 m vertical ground displacement, 14.70 m effective bed change, and 10.90 m recomputed initial water-surface elevation.

[Kajiura (1963)](https://doi.org/10.15083/0000033711) · [Abbate et al. (2024), variable-depth generation](https://doi.org/10.5194/nhess-24-2773-2024)
