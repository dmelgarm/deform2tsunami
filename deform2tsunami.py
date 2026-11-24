"""Sample triangular deformation fields onto WGS84 GeoTIFF pixel centers."""

import argparse
from functools import lru_cache
import math
import os
from pathlib import Path
from tempfile import TemporaryDirectory

import h5py
import matplotlib.tri as mtri
import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.transform import from_origin
from rasterio.warp import transform, transform_bounds


FIELDS = ("vertical", "bed-change", "eta0")
WGS84 = CRS.from_epsg(4326)
STEP = 500       # projected grid spacing, metres
DEPTH_STEP = 250 # reference water-depth interval, metres
MARGIN = 50000   # sampled source margin on each side of a filter tile, metres


class Surface:
    """Original triangle connectivity and the selected static field, in projected metres."""

    def __init__(self, path, field, source_crs):
        with h5py.File(path, "r") as h:
            definition = h.attrs.get("crs_wkt")
            if definition is None and source_crs is None:
                raise ValueError("Missing crs_wkt. Supply the projected metre CRS with source_crs.")
            self.crs = CRS.from_user_input(source_crs if definition is None else definition)
            if source_crs is not None and CRS.from_user_input(source_crs) != self.crs:
                raise ValueError("Supplied source_crs disagrees with the file's crs_wkt.")
            if not self.crs.is_projected or not math.isclose(self.crs.linear_units_factor[1], 1):
                raise ValueError("The source must use projected Cartesian metres.")
            if h.attrs.get("units", "m") not in ("m", b"m") or h.attrs.get("connectivity_index_base", 0) != 0:
                raise ValueError("Expected metre coordinates and zero-based triangle indices.")
            vertices = h["geometry/vertices"][:]
            self.cells = h["viz/topology/cells"][:]
            if field == "eta0":
                if "tsunami/eta0" not in h:
                    raise ValueError("Missing /tsunami/eta0. Use the kajiura command to compute it from displacement.")
                self.values = h["tsunami/eta0"][:]
                if self.values.shape != (len(vertices),) or np.isinf(self.values).any():
                    raise ValueError("Expected one finite eta0 or NaN per vertex.")
            else:
                ds = h["vertex_fields/displacement"]
                if ds.shape != (1, len(vertices), 3):
                    raise ValueError("Expected one static (1, number of vertices, 3) displacement field.")
                self.values = ds[0]
                if not np.isfinite(self.values).all():
                    raise ValueError("Displacement contains non-finite values.")

        if (vertices.ndim != 2 or vertices.shape[1] != 3 or not np.isfinite(vertices).all()
                or self.cells.ndim != 2 or self.cells.shape[1] != 3 or self.cells.dtype.kind not in "iu"
                or self.cells.min() < 0 or self.cells.max() >= len(vertices)):
            raise ValueError("Expected finite XYZ vertices and valid three-node surface triangles.")
        # Merge coincident coordinates only for locating triangles. Displacement may jump
        # across a surface-breaking fault, so interpolation still uses the original IDs.
        xy, ids = np.unique(vertices[:, :2], axis=0, return_inverse=True)
        self.find = mtri.Triangulation(*xy.T, ids.ravel()[self.cells]).get_trifinder()
        self.bounds = (*xy.min(axis=0), *xy.max(axis=0))
        self.tri = vertices[self.cells]
        edges = self.tri[:, 1:, :2] - self.tri[:, :1, :2]
        self.inverse = np.linalg.inv(edges.transpose(0, 2, 1))
        rise = self.tri[:, 1:, 2] - self.tri[:, :1, 2]
        self.slope = np.linalg.solve(edges, rise[..., None])[..., 0]

    def sample(self, x, y, field):
        """Return bed elevation and field at 1-D projected coordinates. Outside is NaN."""
        faces = self.find(x, y)
        inside = faces >= 0
        hit = faces[inside]
        nodes = self.cells[hit]
        offset = np.column_stack((x[inside], y[inside])) - self.tri[hit, 0, :2]
        ab = np.einsum("nij,nj->ni", self.inverse[hit], offset)
        w = np.column_stack((1 - ab.sum(axis=1), ab))
        bed = np.full(len(x), np.nan)
        bed[inside] = self.tri[hit, 0, 2] + (offset * self.slope[hit]).sum(axis=1)
        values = np.full(len(x), np.nan)
        if field == "eta0":
            z = (w * self.values[nodes]).sum(axis=1)
            z[(self.tri[hit, :, 2] >= 0).any(axis=1)] = np.nan
        else:
            u = np.einsum("ni,nij->nj", w, self.values[nodes])
            z = u[:, 2]
            if field == "bed-change":
                z = z - (u[:, :2] * self.slope[hit]).sum(axis=1)
        values[inside] = z
        return bed, values


def kajiura_samples(surface):
    """Return a sampler for the source-depth-banded, rapid-generation approximation."""
    tile_size = 128000  # metres; fixed projected tiles bound memory independently of output resolution
    halo = MARGIN // STEP

    @lru_cache(maxsize=8)
    def filtered_tile(tx, ty):
        # Include both tile-edge nodes for bilinear sampling. Neighbors use the same grid origin.
        ticks = np.arange(-halo, tile_size // STEP + halo + 1) * STEP
        gx, gy = np.meshgrid(tx * tile_size + ticks, ty * tile_size + ticks)
        b, change = surface.sample(gx.ravel(), gy.ravel(), "bed-change")
        wet = b.reshape(gx.shape) < 0
        depth = np.where(wet, -b.reshape(gx.shape), 0)
        source = np.where(wet, change.reshape(gx.shape), 0)
        ny, nx = source.shape
        k = 2*np.pi * np.hypot(np.fft.fftfreq(ny, STEP)[:, None], np.fft.rfftfreq(nx, STEP))
        spectrum = np.zeros(k.shape, dtype=complex)
        for j in range(math.ceil(depth.max() / DEPTH_STEP) + 1):
            # Weight source cells before transforming, using linear interpolation in depth.
            weight = np.maximum(1 - np.abs(depth / DEPTH_STEP - j), 0)
            decay = np.exp(-k * j * DEPTH_STEP)
            spectrum += np.fft.rfft2(source * weight) * (2 * decay / (1 + decay**2))
        eta = np.fft.irfft2(spectrum, s=(ny, nx), axes=(0, 1))
        return eta[halo:-halo, halo:-halo].copy()

    def sample(x, y):
        b, _ = surface.sample(x, y, "vertical")
        wet = np.flatnonzero(b < 0)
        result = np.full(x.shape, np.nan)
        tiles = np.floor_divide(np.column_stack((x[wet], y[wet])), tile_size).astype(int)
        for tx, ty in np.unique(tiles, axis=0):
            pts = wet[(tiles[:, 0] == tx) & (tiles[:, 1] == ty)]
            z = filtered_tile(tx, ty)
            fx, fy = (x[pts] - tx*tile_size)/STEP, (y[pts] - ty*tile_size)/STEP
            ix = np.minimum(fx.astype(int), z.shape[1] - 2)
            iy = np.minimum(fy.astype(int), z.shape[0] - 2)
            wx, wy = fx - ix, fy - iy
            result[pts] = ((1-wy) * ((1-wx)*z[iy, ix] + wx*z[iy, ix+1])
                           + wy * ((1-wx)*z[iy+1, ix] + wx*z[iy+1, ix+1]))
        return result
    return sample


def export_raster(source, output, field, *, arcsec=30, bounds=None, source_crs=None, overwrite=False):
    """Write a float32 geographic GeoTIFF. Bounds are (west, south, east, north) in degrees."""
    if field not in (*FIELDS, "eta0-recomputed"):
        raise ValueError(f"Unknown field: {field}")
    source, output = Path(source), Path(output)
    if output.suffix.lower() not in (".tif", ".tiff") or source.resolve() == output.resolve():
        raise ValueError("Choose a separate .tif or .tiff output file.")
    if output.exists() and not overwrite:
        raise FileExistsError(output)
    if not np.isfinite(arcsec) or arcsec <= 0:
        raise ValueError("arcsec must be finite and positive.")
    surface = Surface(source, field, source_crs)
    if bounds is None:
        bounds = transform_bounds(surface.crs, WGS84, *surface.bounds, densify_pts=51)
    west, south, east, north = bounds
    if not (-180 <= west < east <= 180 and -90 <= south < north <= 90):
        raise ValueError("Expected west < east and south < north within WGS84 lon/lat limits.")
    step = arcsec/3600
    # Anchor pixel edges to multiples of the spacing shared by all exports.
    left, right = math.floor(west/step + 1e-10), math.ceil(east/step - 1e-10)
    bottom, top = math.floor(south/step + 1e-10), math.ceil(north/step - 1e-10)
    width, height = right-left, top-bottom
    affine = from_origin(left*step, top*step, step, step)
    sampler = kajiura_samples(surface) if field == "eta0-recomputed" else None

    # Publish only a complete raster. A failed export leaves an existing output intact.
    with TemporaryDirectory(dir=output.parent) as work:
        pending = Path(work) / "raster.tif"
        with rasterio.open(pending, "w", driver="GTiff", width=width, height=height, count=1,
                           dtype="float32", crs=WGS84, transform=affine, nodata=np.nan,
                           tiled=True, blockxsize=256, blockysize=256, compress="deflate",
                           predictor=3, BIGTIFF="IF_SAFER") as dst:
            dst.set_band_description(1, field)
            dst.set_band_unit(1, "m")
            dst.update_tags(field=field, positive="up", sampling="pixel center", source=source.name)
            if sampler is not None:
                dst.update_tags(working_grid_m=STEP, depth_interval_m=DEPTH_STEP, source_margin_m=MARGIN)
            for _, win in dst.block_windows(1):
                lon = (left + win.col_off + np.arange(win.width) + .5) * step
                lat = (top - win.row_off - np.arange(win.height) - .5) * step
                glon, glat = np.meshgrid(lon, lat)
                x, y = transform(WGS84, surface.crs, glon.ravel(), glat.ravel())
                x, y = np.asarray(x), np.asarray(y)
                if sampler is None:
                    _, data = surface.sample(x, y, field)
                else:
                    data = sampler(x, y)
                dst.write(data.reshape(glon.shape).astype("float32"), 1, window=win)
        if overwrite:
            os.replace(pending, output)
        else:
            os.link(pending, output)  # fails safely if another process just wrote this path
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    raster = commands.add_parser("raster", help="export a deformation or archived sea-surface field")
    raster.add_argument("--field", choices=FIELDS, required=True)
    kajiura = commands.add_parser("kajiura", help="recompute sea-surface elevation from displacement")
    for command in (raster, kajiura):
        command.add_argument("source", type=Path, help="triangular surface HDF5")
        command.add_argument("--out", type=Path, required=True, help="output GeoTIFF")
        command.add_argument("--arcsec", type=float, default=30, help="pixel spacing in arcseconds (default 30)")
        command.add_argument("--bounds", type=float, nargs=4, metavar=("WEST", "SOUTH", "EAST", "NORTH"))
        command.add_argument("--source-crs", help="projected metre CRS for files without crs_wkt")
        command.add_argument("--overwrite", action="store_true", help="replace an existing GeoTIFF")
    args = parser.parse_args()
    field = "eta0-recomputed" if args.command == "kajiura" else args.field
    print(export_raster(args.source, args.out, field, arcsec=args.arcsec, bounds=args.bounds,
                        source_crs=args.source_crs, overwrite=args.overwrite))


if __name__ == "__main__":
    main()
