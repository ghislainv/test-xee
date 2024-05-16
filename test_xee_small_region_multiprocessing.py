"""Test for using xee."""

import os
import multiprocessing as mp
from glob import glob

import ee
import xarray as xr
import forestatrisk as far
from osgeo import gdal

from xarray2geotiff_2 import save_xarray
from make_grid import create_buffer, make_grid, grid_intersection

iso = "COD"
years = [2000, 2005, 2010, 2015, 2020]

# Initialize Earth Engine
ee.Initialize(project="forestatrisk",
              opt_url=("https://earthengine-highvolume"
                       ".googleapis.com"))

# Download borders
far.data.download.download_gadm(iso, output_dir="data")
borders_shp = os.path.join("data", f"gadm36_{iso}_0.shp")

# Buffer of 10km around borders
buff = 0.08983152841195216  # in dd, ~10 km
buff_file = os.path.join("data", f"gadm36_{iso}_buffer.gpkg")
create_buffer(input_file=borders_shp,
              output_file=buff_file,
              buffer_dist=buff)

# Compute extent
extent_latlong = far.get_vector_extent(buff_file)

# region = ee.Geometry.Rectangle(extent_latlong,
#                                proj=proj,
#                                geodesic=False)
# region = region.buffer(10000).bounds()

# # Grid
# res = 0.5  # in degrees
# projection = ee.Projection(proj).scale(res, res)
# grid = region.coveringGrid(projection)
# nfeatures = grid.size().int()
# new_region = grid.geometry().bounds()

# ========================
# Make minimal grid
# ========================

proj = "EPSG:4326"
epsg_code = 4326
csize = 1.0  # 1 decimal degree (dd)
scale = 0.000269494585235856472  # in dd, ~30 m
grid_gpkg = os.path.join("outputs", f"grid_{iso}.gpkg")
grid = make_grid(
    extent_latlong,
    buff=0,
    csize=csize,
    scale=scale,
    proj=epsg_code,
    ofile=grid_gpkg)
ofile = os.path.join("outputs", f"mingrid_{iso}.gpkg")
grid_i = grid_intersection(grid, grid_gpkg, ofile, borders_shp)

# =========================
# Forest cover change data
# =========================

AP = ee.ImageCollection("projects/JRC/TMF/v1_2022/AnnualChanges")
AP = AP.mosaic().toByte()

# ap_all_year: forest if Y = 1 or 2.
ap_forest = AP.where(AP.eq(2), 1)
ap_all_year = ap_forest.where(ap_forest.neq(1), 0)

forest_list = []
band_names = []

for year in years:
    id_year = year - 1990 - 1
    ap = ap_all_year.select(list(range(id_year, 33)))
    forest_yr = ap.reduce(ee.Reducer.sum()).gte(1)
    forest_yr = forest_yr.set(
        "system:time_start",
        ee.Date.fromYMD(year, 1, 1).millis())
    forest_list.append(forest_yr)
    band_names.append(f"forest{year}")

forest_ic = ee.ImageCollection(forest_list)


def get_date(image):
    """Get formatted date."""
    date = ee.Image(image).date().format("YYYY-MM-dd")
    return date


def filter_and_mosaic(d):
    """Create mosaic for one date."""
    d = ee.Date(d)
    im = (forest_ic
          .filterDate(d, d.advance(1, "day"))
          .mosaic().toByte())
    im = im.set("system:time_start", d.millis(),
                "system:id", d.format("YYYY-MM-dd"))
    return im


def mosaic_by_date(img_list):
    """Mosaic by date."""
    unique_dates = img_list.map(get_date).distinct()
    mosaic_list = unique_dates.map(filter_and_mosaic)
    return ee.ImageCollection(mosaic_list)


forest = mosaic_by_date(ee.List(forest_list))


# ds = (
#     xr.open_dataset(
#         forest,
#         engine="ee",
#         crs=proj,
#         scale=scale,
#         chunks={},
#         geometry=grid_i[0],
#     )
#     .astype("b")
#     .rename({"lon": "longitude", "lat": "latitude"})
#     .rename({"sum": "forest_cover"})
# )


def write_geotiff(index, extent):
    """Write geotiff."""

    # Open dataset
    ds = (
        xr.open_dataset(
            forest,
            engine="ee",
            crs=proj,
            scale=scale,
            geometry=extent,
        )
        .astype("b")
        .rename({"lon": "longitude", "lat": "latitude"})
        .rename({"sum": "forest_cover"})
    )

    # Load and write data to geotiff
    fname = f"outputs/forest_{iso}.tif"
    save_xarray(ds, "forest_cover", fname, index)


# Multiprocessing
import time
start_time = time.time()
ncpu = os.cpu_count() - 1
pool = mp.Pool(processes=ncpu)
args = [(i, ext) for (i, ext) in enumerate(grid_i)]
pool.starmap_async(write_geotiff, args).get()
pool.close()
pool.join()
stop_time = time.time()
print("--- %s seconds ---" % (stop_time - start_time))

# Make vrt
tif_forest_files = glob("outputs/forest_" + iso + "*.tif")
# Callback
verbose = True
cback = gdal.TermProgress if verbose else 0
vrt_dataset = gdal.BuildVRT("outputs/forest.vrt",
              tif_forest_files,
              callback=cback)
vrt_dataset.FlushCache()

# End Of File
