import fsspec
import xarray as xr
url = "https://gamma.hdfgroup.org/ftp/pub/outgoing/NASAHDF/ATL06_20190223232535_08780212_001_01.h5"
with fsspec.open(url) as f:
    ds = xr.open_dataset(f)
    print(ds)