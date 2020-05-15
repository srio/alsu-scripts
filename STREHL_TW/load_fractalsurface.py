import numpy

a = numpy.loadtxt("fractalsurface.dat")

print(a.shape)

from srxraylib.plot.gol import plot_image, plot

x = a[1:,0].copy()
y = a[0,1:].copy()
z = a[1:,1:].copy()

print(z.shape,x.shape,y.shape)
plot_image(z, x, y,  aspect='auto')
zcentral = (z[:, z.shape[1] // 2]).copy()

plot(1e3 * x, 1e9 * zcentral )
print( zcentral.std())

f = open("fractalprofile.dat",'w')
for i in range(zcentral.size):
    f.write("%g  %g\n" % (x[i], zcentral[i]))
f.close()
print("File written to disk: fractalprofile.dat")


# print(a[0,1:], a[1:,0],)