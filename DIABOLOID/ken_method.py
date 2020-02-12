


import numpy

p = 29.3
q = 19.53
theta = 4.5e-3

y = numpy.linspace(-0.1, 0.1, 100)
x = numpy.linspace(-0.01,0.01, 80)

X = numpy.outer(x,numpy.ones_like(y))
Y = numpy.outer(numpy.ones_like(x),y)


s = p * numpy.cos(2*theta)
z0 = p * numpy.sin(2*theta)
c = p + q

Z = - numpy.sqrt( c**2 + q**2 - s**2 - 2 * Y * (s + q) - 2 * c * numpy.sqrt( X**2 + (q - Y) ) )

print(Z.shape)
from srxraylib.plot.gol import plot_image

plot_image(Z) #,X,Y)
