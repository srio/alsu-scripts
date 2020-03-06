

import numpy
# from oasys.util.oasys_objects import OasysSurfaceData
from oasys.util.oasys_util import write_surface_file



def get_radius(method="howard"):
    if method == "howard":
        # howard
        h = 15e-6
        x = 10e-3
    elif method == "X":
        # Along X
        h = 3e-6
        x = 15e-3
    elif method == "Y":
        # Along Y
        h = 5e-6
        x = 10e-3

    R = (h ** 2 + x ** 2) / (2 * h)
    gamma = 2 * h * x / (x ** 2 + h ** 2)
    print("method: %s , Radius: %f, gamma: %f "%(method,R,gamma))

    return R


Ly = 30e-3
Lx = 20e-3

y = numpy.linspace(-Ly/2,Ly/2,201)
x = numpy.linspace(-Lx/2,Lx/2,41)

radius_y = -get_radius("X")
radius_x =  get_radius("Y")

# radius_y = -get_radius("howard")
# radius_x =  get_radius("howard") #/ 10

# radius_y = -1000.0
# radius_x = 1000.0

if radius_x > 0:
    z_x = radius_x - numpy.sqrt(radius_x ** 2 - x ** 2)
else:
    z_x = radius_x + numpy.sqrt(radius_x ** 2 - x ** 2)

# height = radius + numpy.sqrt(radius ** 2 - x2_oe ** 2)
if radius_y < 0:
    z_y = radius_y + numpy.sqrt(radius_y ** 2 - y ** 2)
else:
    z_y = radius_y - numpy.sqrt(radius_y ** 2 - y ** 2)

from srxraylib.plot.gol import plot, plot_surface

# z_x *= 0
# z_y *= 0
plot(1e3*x,1e9*z_x,1e3*y,1e9*z_y,legend=["z vs x","z vs y"],xtitle="x or y coordinate [mm]",ytitle="Height [nm]")

# Z = numpy.outer(z_x,z_y)
Z = numpy.zeros( (x.size, y.size) )

for i in range(x.size):
    for j in range(y.size):
        Z[i,j] = z_x[i] + z_y[j]

plot_surface(1e9*Z,1e3*x,1e3*y,xtitle="Shadow X [mm]",ytitle="Shadow Y [mm]",ztitle="Height Z [nm]",title="")

# OasysSurfaceData(xx=x,
#                  yy=y,
#                  zz=Z,
#                  surface_data_file="BL822_crystal_error.h5")

write_surface_file(Z.T, x, y, "BL822_crystal_error.h5", overwrite=True)
print("write_h5_surface: File for OASYS " + "BL822_crystal_error.h5" + " written to disk.")


