

import h5py, numpy
from srxraylib.plot.gol import plot_image, plot_scatter, plot
from srxraylib.util.inverse_method_sampler import Sampler2D
import Shadow

f = h5py.File("manolone.hdf5",'r')


Z = f["/surface_file/Z"][:]
x = f["/surface_file/X"][:]
y = f["/surface_file/Y"][:]
f.close()

# plot_image(Z,x,y)

s = Sampler2D(numpy.rot90(Z),x,y)
s._cdf_calculate()

print(s._cdf1.shape,s._cdf2.shape)

# plot_image(s._cdf2)


npoints = 500000
points_x, points_y = s.get_n_sampled_points(npoints)

# print(points.shape)

# plot_scatter(points_x,points_y)

oe0 = Shadow.Source()

oe0.NPOINT = npoints

beam = Shadow.Beam()

beam.genSource(oe0)
beam.rays[:,0] = points_x
beam.rays[:,2] = points_y

beam.write("begin.dat")

beam.histo2(1,3,)


#
# nruns = 100
#
# xm = numpy.zeros(nruns)
# ym = numpy.zeros(nruns)
# im = numpy.arange(nruns)
#
# for i in range(nruns):
#     s = Sampler2D(Z, x, y)
#     points_x, points_y = s.get_n_sampled_points(100000)
#     xm[i] = numpy.average(points_x)
#     ym[i] = numpy.average(points_y)
#
# plot(im,xm,im,ym)
#
# print("X: %g +/-  %g"%(numpy.average(xm),numpy.std(xm)))
# print("Y: %g +/-  %g"%(numpy.average(ym),numpy.std(ym)))








