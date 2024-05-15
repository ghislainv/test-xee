"""Test for using xee."""

import ee
import numpy as np
import xarray as xr
import forestatrisk as far

from xarray2geotiff_2 import save_xarray

iso3 = "REU"
years = [2000, 2005, 2010, 2015, 2020]

ee.Initialize(project="forestatrisk",
              opt_url=("https://earthengine-highvolume"
                       ".googleapis.com"))

far.data.download.download_gadm(iso3, output_dir="data")
extent_latlong = far.extent_shp(f"data/gadm36_{iso3}_1.shp")

# Region
proj = "EPSG:4326"
region = ee.Geometry.Rectangle(extent_latlong,
                               proj=proj,
                               geodesic=False)
region = region.buffer(10000).bounds()

# Grid
res = 0.5  # in degrees
projection = ee.Projection(proj).scale(res, res)
grid = region.coveringGrid(projection)
new_region = grid.geometry().bounds()

# Scale (in decimal degree)
scale = 0.0002727272727

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

# forest = forest.toBands()

# forest_clip = forest.clipToBoundsAndScale(
#     geometry=grid.first().geometry(),
#     scale=30)

# Open dataset
ds = (
    xr.open_dataset(
        forest,
        engine="ee",
        crs=proj,
        scale=scale,
        chunks={"time": -1, "lat": 1024, "lon": 1024},
        geometry=region,
    )
    .astype("b")
    .rename({"lon": "longitude", "lat": "latitude"})
    .rename({"sum": "forest_cover"})
)


# # This is working !
# xarr = ds.load()
# data_var = "forest_cover"
# fname = f"outputs/forest_{iso3}.tif"
# save_xarray(ds, "forest_cover", fname=fname)


# # Open dataset
# ds = (
#     xr.open_dataset(
#         forest,
#         engine="ee",
#         crs=proj,
#         scale=scale,
#         chunks={"time": -1, "lat": 512, "lon": 512},
#         geometry=grid.first().geometry(),
#     )
#     .astype("b")
#     .rename({"lon": "longitude", "lat": "latitude"})
#     .rename({"sum": "forest_cover"})
# )

# Make template for output xarray
nchunks_lat = len(ds.chunksizes["latitude"])
nchunks_lon = len(ds.chunksizes["longitude"])
data_fc = np.zeros([1, nchunks_lat, nchunks_lon])
ds_template = xr.Dataset(
    data_vars={
        "forest_cover": (("time", "latitude", "longitude"), data_fc)
    },
    coords={
        "time": ("time", [0]),
        "longitude": ("longitude", list(range(nchunks_lon))),
        "latitude": ("latitude", list(range(nchunks_lat))),
    },
)
# Set chunk size (not number of chunks !!!)
ds_template = (ds_template
               .chunk(chunks={"time": 1,
                              "latitude": 1,
                              "longitude": 1}))

# ds_template = xr.ones_like(ds.time).chunk({"time": 5})

# ds_template = (
#     xr.ones_like(ds[["time", "latitude", "longitude"]])
#     .chunk({"latitude": 1024, "longitude": 1024})
# )

nchunks_lat = len(ds.chunksizes["latitude"])
nchunks_lon = len(ds.chunksizes["longitude"])
nchunks = nchunks_lat * nchunks_lon
ds_template = xr.Dataset(
    data_vars={
        "value": (("block"), [0] * nchunks)
    },
    coords={
        "block": ("block", list(range(nchunks))),
    },
).chunk({"block": 1}).value

# Using map_blocks
ofile = f"outputs/forest_{iso3}.tif"
tasks = ds.map_blocks(
    save_xarray,
    kwargs={"data_var": "forest_cover",
            "fname": ofile},
    template=ds_template).compute()

# tasks.compute(num_workers=4, scheduler="threads")

# End Of File
