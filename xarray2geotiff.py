"""Saves an xarray dataset to a Cloud Optimized GeoTIFF (COG).
https://gist.github.com/GerardoLopez/35123d4a15aa31f3ea4b01efb5b26d4d
"""


# import ee
import os
from osgeo import gdal, gdal_array
from osgeo import osr


# ee.Initialize(project="forestatrisk",
#               opt_url=("https://earthengine-highvolume"
#                        ".googleapis.com"))


def get_dst_dataset(dst_img, cols, rows, layers, dtype, proj, gt):
    """
    Create a GDAL data set in Cloud Optimized GeoTIFF (COG) format
    :param dst_img: Output filenane full path
    :param cols: Number of columns
    :param rows: Number of rows
    :param layers: Number of layers
    :param dtype: GDAL type code
    :param proj: Projection information in WKT format
    :param gt: GeoTransform tupple
    :return dst_ds: GDAL destination dataset object
    """
    gdal.UseExceptions()
    try:
        # Default driver options to create a COG
        driver = gdal.GetDriverByName('GTiff')
        driver_options = ['COMPRESS=DEFLATE',
                          'PREDICTOR=1',
                          'BIGTIFF=YES',
                          'TILED=YES',
                          'COPY_SRC_OVERVIEWS=YES']

        # Create driver
        dst_ds = driver.Create(dst_img, cols, rows, layers,
                               dtype, driver_options)

        # Set cartographic projection
        dst_ds.SetProjection(proj)
        dst_ds.SetGeoTransform(gt)

    except Exception as err:
        if err.err_level >= gdal.CE_Warning:
            # print('Cannot write dataset: %s' % self.input.value)
            # Stop using GDAL exceptions
            gdal.DontUseExceptions()
            raise RuntimeError(err.err_level, err.err_no, err.err_msg)

    gdal.DontUseExceptions()
    return dst_ds


def get_resolution_from_xarray(xarray):
    """
    Method to create a tuple (x resolution, y resolution) in x and y
    dimensions.

    :param xarray: xarray with latitude and longitude variables.

    :return: tuple with x and y resolutions
    """

    x_res = xarray.longitude.values[1] - xarray.longitude.values[0]
    y_res = xarray.latitude.values[1] - xarray.latitude.values[0]

    return (x_res, y_res)


def save_xarray(xarray, data_var, fname):
    """
    Saves an xarray dataset to a Cloud Optimized GeoTIFF (COG)
    :param xarray: xarray Dataset
    :param data_var: Data variable in the xarray dataset
    :param fname: Full path of file where to save the data
    """

    # Create GeoTransform - perhaps the user requested a
    # spatial subset, therefore is compulsory to update it

    # GeoTransform -- case of a "north up" image without
    #                 any rotation or shearing
    #  GeoTransform[0] top left x
    #  GeoTransform[1] w-e pixel resolution
    #  GeoTransform[2] 0
    #  GeoTransform[3] top left y
    #  GeoTransform[4] 0
    #  GeoTransform[5] n-s pixel resolution (negative value)

    # Create tmp xarray DataArray
    _xarray = getattr(xarray, data_var)
    print(_xarray)

    x_res, y_res = get_resolution_from_xarray(_xarray)

    gt = (_xarray.longitude.data[0] - (x_res / 2.),
          x_res, 0.0,
          _xarray.latitude.data[0] - (y_res / 2.),
          0.0, y_res)

    # Coordinate Reference System (CRS) in a PROJ4 string to a
    # Spatial Reference System Well Known Text (WKT)
    crs = xarray.attrs['crs']
    crs = int(crs[5:9])
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(crs)
    proj = srs.ExportToWkt()

    # Get GDAL datatype from NumPy datatype
    if _xarray.dtype == 'bool':
        dtype = gdal.GDT_Byte
    else:
        dtype = gdal_array.NumericTypeCodeToGDALTypeCode(_xarray.dtype)

    # Dimensions
    layers, rows, cols = _xarray.shape

    # File name with postfix
    postfix = f"{int(gt[0])}_{int(gt[3])}"
    dirname = os.path.dirname(fname)
    fname = os.path.splitext(os.path.basename(fname))[0]
    fname = os.path.join(dirname, f"{fname}_{postfix}.tif")

    # Create destination dataset
    dst_ds = get_dst_dataset(
        dst_img=fname, cols=cols, rows=rows,
        layers=layers, dtype=dtype, proj=proj, gt=gt)

    for layer in range(layers):
        dst_band = dst_ds.GetRasterBand(layer + 1)

        # Date
        if 'time' in _xarray.dims:
            dst_band.SetMetadataItem(
                'time',
                _xarray.time.data[layer].astype(str))

        # Data variable name
        dst_band.SetMetadataItem('data_var', data_var)

        # Data
        dst_band.WriteArray(_xarray[layer].data)

    dst_band = None
    del dst_ds

    return xarray

# End Of File
