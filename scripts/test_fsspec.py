import fsspec
import xarray as xr
import h5py
url = "https://github.com/srio/alsu-scripts/blob/master/scripts/ISN_KBH.hdf5?raw=true"
with fsspec.open(url) as f:
    ds = xr.open_dataset(f)
    print(ds,type(ds))
    # z = ds["surface_file/Z"]
    # print(z.shape)



# import numpy
# import json
# arr = numpy.linspace(-10, 10, 100).reshape((10,10))
# print(json.dumps({'mydata': arr.tolist()}))