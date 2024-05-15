"""Test for using xee."""

import ee
import xarray as xr
import forestatrisk as far

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
res = 1  # in degrees
projection = ee.Projection(proj).scale(res, res)
grid = region.coveringGrid(projection)
new_region = grid.geometry().bounds()

# Scale (in decimal degree)
scale = 0.0002727272727

AP = ee.ImageCollection("projects/JRC/TMF/v1_2022/AnnualChanges")
AP = AP.mosaic().toByte()


# JRC annual product (AP)
def compute_forest(index):
    """Compute forest."""
    # Select feature
    feat_coll = grid.filter(ee.Filter.eq('system:index', index))
    ap_clip = AP.clip(feat_coll)

    # ap_all_year: forest if Y = 1 or 2.
    ap_clip_forest = ap_clip.where(ap_clip.eq(2), 1)
    ap_all_year = ap_clip_forest.where(ap_clip_forest.neq(1), 0)

    # forest = ee.Image()
    _forest_list = []
    band_names = []

    for year in years:
        id_year = year - 1990 - 1
        ap = ap_all_year.select(list(range(id_year, 33)))
        forest_yr = ap.reduce(ee.Reducer.sum()).gte(1)
        forest_yr = forest_yr.set(
            "system:time_start",
            ee.Date.fromYMD(year, 1, 1).millis())
        # forest = forest.addBands(forest_yr)
        _forest_list.append(forest_yr)
        band_names.append(f"forest{year}")

    # forest = forest.select([i for (i, _) in enumerate(years)], band_names)
    # forest = forest.set("system:bandNames", band_names)
    # forest = forest

    return _forest_list


# Make computations
index_list = grid.aggregate_array('system:index')
forest_list = index_list.map(compute_forest).flatten()
# forest = ee.ImageCollection.fromImages(forest_list).mosaic()
# forest = ee.ImageCollection(forest)
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


forest = mosaic_by_date(forest_list)

# Open dataset
ds = (
    xr.open_dataset(
        forest,
        engine="ee",
        crs=proj,
        scale=scale,
        chunks={},
        geometry=new_region
    )
    .astype("b")
    .rename({"time": "bands", "lon": "x", "lat": "y"})
    .rename({"sum": "forest_cover"})
    .transpose("bands", "y", "x")
)

forest_rio = (
    ds["forest_cover"]
    .rio.set_spatial_dims(x_dim="x", y_dim="y")
    .rio.write_crs("epsg:4326")
    .rio.write_coordinate_system()
)

# Save to disk
forest_rio.rio.to_raster(
    f"outputs/forest_{iso3}.tif",
    tiled=True,
    driver="GTiff",
    compress="LZW"
)

# End Of File
