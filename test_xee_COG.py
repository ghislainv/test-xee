"""Test for using xee."""

import ee
import xarray as xr
import pandas as pd
import rioxarray
from geefcc.download_gadm import download_gadm
from geefcc.get_vector_extent import get_vector_extent
from dask.distributed import LocalCluster, Client

#from xarray2geotiff import save_xarray

iso3 = "PER"
years = [2000, 2005, 2010, 2015, 2020]

ee.Initialize(project="forestatrisk",
              opt_url="https://earthengine-highvolume.googleapis.com")

ofile = f"data/gadm41_{iso3}_0.gpkg"
download_gadm(iso3, output_file=ofile)
extent_latlong = get_vector_extent(ofile)

# Region
region = ee.Geometry.Rectangle(extent_latlong,
                               proj="EPSG:4326",
                               geodesic=False)
region = region.buffer(10000).bounds()

# Scale
scale = 0.0002727272727

# JRC annual product (AP)
IC = ee.ImageCollection("projects/JRC/TMF/v1_2022/AnnualChanges")
IC = IC.mosaic().toByte()
IC = IC.select(list(range(9, 33)))

# Open dataset, IC should be an image collection.
IC = ee.ImageCollection(IC)
ds = (xr.open_dataset(
    IC,
    engine="ee",
    crs="EPSG:4326",
    scale=scale,
    chunks={},
    geometry=region)
      .astype("b")
      .drop_vars("time")
      .squeeze())

# Reorganize xarray.Dataset into xarray.DataArray
dec_years_list = list(range(years[0]-1, 2023))
years_list = [i + 1 for i in dec_years_list]
AP = (xr.concat(
    [ds[f"Dec{i}"] for i in dec_years_list],
    dim=pd.Index(years_list, name="years"))
      .rename("Annuals Products")
      .chunk({"years": -1, "lat": 256, "lon": 256}))


# AP_allYear: forest if Y = 1 or 2.
def forest_masking(annual_product):
    """Masking forest."""
    ap_forest = xr.where((annual_product == 1)
                         + (annual_product == 2), 1, 0)
    return ap_forest


# Forest cover for each year
def determine_forest(ap_forest, year, year_end):
    """Determine forest."""
    return (
        ap_forest
        .sel(years=slice(year, year_end + 1))
        .sum("years")
        .rename("forest" + str(year)) >= 1
    )


AP_forest = AP.map_blocks(forest_masking)
forest = xr.merge(
    [AP_forest.map_blocks(
        determine_forest,
        kwargs={"year": year, "year_end": 2022}
    ) for year in years]
)

# Reorganize
proj = "EPSG:4326"
forest_years_list = list(range(2000, 2025, 5))
forest_rio = xr.concat(
    [forest[f"forest{i}"] for i in forest_years_list],
    dim=pd.Index(forest_years_list, name="bands"))
forest_rio = forest_rio.chunk({"bands": -1, "latitude": 512, "longitude": 512})

forest_rio2 = (forest_rio
              .astype('b')
              .rename({"longitude": "x", "latitude": "y"})
              .transpose("bands", "y", "x")
              .rename("forest"))

# Export to GeoTiff
# save_xarray(fname=f"outputs/fcc_{iso3}.tif",
#             xarray=ds, data_var="fcover")

# Save to disk
with LocalCluster() as cluster, Client(cluster) as client:
    forest_rio2.rio.to_raster(
        f"outputs/forest_{iso3}.tif",
        tiled=True,
        driver="GTiff",
        compress="LZW",
        # lock=Lock("rio", client=client),
    )

# End Of File
