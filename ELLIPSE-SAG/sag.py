import numpy
from shadow4.optical_surfaces.conic import Conic
from srxraylib.plot.gol import plot

################################
# Inputs
################################
# p = 15
# q = 5
# theta1 = 2e-3
# mirror_length = 0.4
number_of_points = 300

theta1 = 21.817e-3
mirror_length = 175e-3
p = 28.6
q = 5.23
###################################

coeffs_ellipse = Conic.initialize_as_ellipsoid_from_focal_distances(p, q, theta1, cylindrical=1, cylangle=0,)

coeffs_cylinder = Conic.initialize_as_sphere_from_focal_distances(p, q, theta1, cylindrical=1, cylangle=0)

print( coeffs_ellipse.info() )

# OE surface in form of conic equation:
#   ccc[0]*X^2 + ccc[1]*Y^2 + ccc[2]*Z^2
#   ccc[3]*X*Y + ccc[4]*Y*Z + ccc[5]*X*Z
#   ccc[6]*X   + ccc[7]*Y   + ccc[8]*Z + ccc[9] = 0
# we need Z vs Y (X=0):
#   ccc[2]*Z^2
#   + ccc[4]*Y*Z + ccc[8]*Z
#   + ccc[7]*Y   ccc[1]*Y^2 + ccc[9] = 0
#
#  Z is given by second degree equation solved in the method height:

y = numpy.linspace(-0.5 * mirror_length, 0.5 * mirror_length, number_of_points)

z_ellipse = coeffs_ellipse.height(y=y,x=0,return_solution=0)
z_cylinder = coeffs_cylinder.height(y=y,x=0,return_solution=0)


plot(y, 1e6 * z_ellipse,
     y, 1e6 * z_cylinder,
     legend=["ellipse", "cylinder"],
     xtitle="y [m]", ytitle="Z [um]",
     title="p=%5.3f m, q=%5.3f m, theta=%5.3f mrad" % (p, q, 1e3 * theta1))

n0 = number_of_points * 3 // 4
plot(y, 1e9 * (z_ellipse - z_cylinder),
     legend=["ellipse - cylinder"],
     xtitle="y [m]", ytitle="Zdiff [nm]",
     title="p=%5.3f m,  q=%5.3f m, theta=%5.3f mrad" % (p, q, 1e3 * theta1))



coeffs = numpy.polyfit(y[n0:], 1e9 * (z_ellipse - z_cylinder)[n0:], 2, rcond=None, full=False, w=None, cov=False)
yfit = coeffs[2] + y * coeffs[1] + y**2 * coeffs[0]
# yfit = coeffs[1] + y * coeffs[0]
print("Coeffs: ", coeffs)

plot(y, 1e9 * (z_ellipse - z_cylinder),
     y[n0:], (yfit)[n0:],
     legend=["ellipse - cylinder", "quadratic fit"],
     xtitle="y [m]", ytitle="Zdiff [nm]",
     title="p=%5.3f m,  q=%5.3f m, theta=%5.3f mrad" % (p, q, 1e3 * theta1))