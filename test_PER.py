"""Test for using xee."""

import time

import ee
from geefcc.download_gadm import download_gadm
from geefcc.get_vector_extent import get_vector_extent

from xee_jrc import run_task

ee.Initialize(project="forestatrisk",
              opt_url="https://earthengine-highvolume.googleapis.com")

iso3 = "PER"
years = [2000, 2010, 2020]

border_file = f"data/gadm41_{iso3}_0.gpkg"
download_gadm(iso3, output_file=border_file)
extent_latlong = get_vector_extent(border_file)
scale = 0.000269494585235856472  # in dd, ~30 m
proj = "EPSG:4326"

start_time = time.time()
run_task(iso3, years, extent_latlong, scale, proj)
end_time = time.time()

ellapsed_time = (end_time - start_time) / 60
start_time = time.time()

# ########
# """Using GEE to get forest cover change from TMF."""

# # Annual product legend
# # ---------------------
# # 1. Undisturbed Tropical moist forest (TMF)
# # 2. Degraded TMF
# # 3. Deforested land
# # 4. Forest regrowth
# # 5. Permanent or seasonal water
# # 6. Other land cover

# # Third party imports
# import ee
# import xarray as xr
# import pandas as pd



# # Region
# region = ee.Geometry.Rectangle(extent_latlong,
#                                proj="EPSG:4326",
#                                geodesic=False)
# region = region.buffer(10000).bounds()

# # JRC annual product (AP)
# IC = ee.ImageCollection("projects/JRC/TMF/v1_2022/AnnualChanges")
# ds = xr.open_dataset(IC,
#                      engine="ee",
#                      crs="EPSG:4326",
#                      scale=scale,
#                      geometry=region,
#                      ).isel(time=2).drop_vars('time').squeeze()
# ds = ds.chunk({"lon": 1500, "lat": 1500})
# ds

# # Reorganize AP dataset
# date_list = range(years[0]-1, years[-1] + 1)
# AP = xr.concat(
#     [ds[f"Dec{i}"] for i in date_list],
#     dim=pd.Index(date_list, name="years")
#     ).rename("Annuals Products")

# # AP_allYear: forest if Y = 1 or 2.
# def forest_masking(AP):
#     return xr.where((AP == 1) + (AP == 2), 1, 0)
# AP_forest = AP.map_blocks(forest_masking)

# # Forest cover for each year
# def determine_forest(AP_forest, year, ref_year):
#     return AP_forest.\
#         sel(years=slice(year-1, ref_year+1)).\
#         sum("years").\
#         rename("forest"+str(year)) >= 1

# forest = xr.merge(
#     [AP_forest.map_blocks(
#         determine_forest,
#         kwargs={"year": year, "ref_year": years[-1]}
#         ) for year in years]
#     )

# # Prepare the raster and save it to disk
# forest_date_list = [f"forest{i}" for i in years]
# forest = xr.concat(
#     [forest[i] for i in forest_date_list],
#     dim=pd.Index(forest_date_list, name="bands"))
# forest = (forest
#           .astype('b')
#           .rio.set_spatial_dims(x_dim="lon", y_dim="lat")
#           .rio.write_crs("epsg:4326")
#           .rio.write_coordinate_system()
#           .rename({"lon": "x", "lat": "y"})
#           .transpose("bands", "y", "x")
#           .rename("forest"))

# # Save to disk
# forest.rio.to_raster(
#     f"forest_{iso3}.tif",
#     tiled=True,
#     driver="GTiff",
#     compress="LZW"
# )


# End Of File
