


import numpy
from oasys.util.oasys_util import write_surface_file


p = 29.3
q = 19.53
theta = 4.5e-3

y = numpy.linspace(-0.1, 0.1, 1000)
x = numpy.linspace(-0.01,0.01, 200)

X = numpy.outer(x,numpy.ones_like(y))
Y = numpy.outer(numpy.ones_like(x),y)


s = p * numpy.cos(2*theta)
z0 = p * numpy.sin(2*theta)
c = p + q

print("s: %f, z0: %f, c: %f"%(s,z0,c))

Z = - numpy.sqrt( c**2 + q**2 - s**2 - 2 * Y * (s + q) - 2 * c * numpy.sqrt( X**2 + (q - Y)**2 ) )

Z += z0 # Z[Z.shape[0]//2,Z.shape[1]//2]
print(Z.shape,Z.min(),Z.max())
from srxraylib.plot.gol import plot_image

plot_image(Z,1e3*x,1e3*y,aspect='auto')

write_surface_file(Z.T, x, y, "C:\\Users\\manuel\\Oasys\\ken_method.h5", overwrite=True)
