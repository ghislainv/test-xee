"""Test for using xee."""

import ee
import xarray as xr
import forestatrisk as far

from xarray2geotiff import save_xarray

iso3 = "MTQ"
years = [2000, 2005, 2010, 2015, 2020]

ee.Initialize(project="forestatrisk",
              opt_url="https://earthengine-highvolume.googleapis.com")

far.data.download.download_gadm(iso3, output_dir="data")
extent_latlong = far.extent_shp("data/gadm36_MTQ_1.shp")

# Region
region = ee.Geometry.Rectangle(extent_latlong,
                               proj="EPSG:4326",
                               geodesic=False)
region = region.buffer(10000).bounds()

# Scale
scale = 0.0002727272727

# JRC annual product (AP)
AP = ee.ImageCollection("projects/JRC/TMF/v1_2020/AnnualChanges")
AP = AP.mosaic().toByte()
AP = AP.clip(region)

# ap_allYear: forest if Y = 1 or 2.
AP_forest = AP.where(AP.eq(2), 1)
ap_allYear = AP_forest.where(AP_forest.neq(1), 0)

# Forest in Jan 2020
ap_2020_2021 = ap_allYear.select(list(range(29, 31)))
forest2020 = ap_2020_2021.reduce(ee.Reducer.sum()).gte(1)

# Forest cover Jan 2015
ap_2015_2021 = ap_allYear.select(list(range(24, 31)))
forest2015 = ap_2015_2021.reduce(ee.Reducer.sum()).gte(1)

# Forest cover Jan 2010
ap_2010_2021 = ap_allYear.select(list(range(19, 31)))
forest2010 = ap_2010_2021.reduce(ee.Reducer.sum()).gte(1)

# Forest cover Jan 2005
ap_2005_2021 = ap_allYear.select(list(range(14, 31)))
forest2005 = ap_2005_2021.reduce(ee.Reducer.sum()).gte(1)

# Forest cover Jan 2000
ap_2000_2021 = ap_allYear.select(list(range(9, 31)))
forest2000 = ap_2000_2021.reduce(ee.Reducer.sum()).gte(1)

# Forest raster with five bands
forest = ee.ImageCollection([forest2000, forest2005, forest2010,
                             forest2015, forest2020])


# Set time_start
def set_time(image):
    """Set time."""
    return image.set(
        'system:time_start',
        ee.Date.fromYMD(2000, 1, 1).millis())


forest_date = forest.map(set_time)

# Open dataset
ds = xr.open_dataset(forest_date,
                     engine="ee",
                     crs="EPSG:4326",
                     scale=scale,
                     chunks={},
                     geometry=region).astype("b")

# Export to GeoTiff
save_xarray(fname=f"outputs/fcc_{iso3}.tif",
            xarray=ds, data_var="fcover")

# End Of File
