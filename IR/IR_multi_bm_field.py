import numpy
from srxraylib.plot.gol import plot, set_qt
from scipy.ndimage import gaussian_filter1d

set_qt()

L = 1605.0 #mm

y = numpy.linspace(0,L, 2000)

B = y * 0.0


for i in range(y.size):
    if y[i] > 75 and y[i] < 575: B[i] = -0.876
    if y[i] > 650 and y[i] < 975: B[i] = 0.16
    if y[i] > 1030 and y[i] < 1530: B[i] = -0.8497

# plot(y, B)


B2 = gaussian_filter1d(B, 2.5)

yy = y.copy()
yy -= yy[y.size//2]
yy *= 1e-3

# plot(yy, B, yy, B2, legend=["original","smoothed"],xtitle="y / m",ytitle="B / T")

f = open("BM_multi.b", "w")
for i in range(y.size):
    f.write("%f  %f\n" % (yy[i], B2[i]))
f.close()
print("File written to disk: BM_multi.b")



# analyse M1

B2 = B.copy()
B2[numpy.where(y > 575)] = 0.0
electron_energy_in_GeV = 1.9
# plot(yy,B2)
radius = 3.334728*electron_energy_in_GeV/B
w = B2 / B2.min()
# plot(yy,w)
center = numpy.average(yy,weights=w)

t = numpy.where(w > 0.5)
width = yy[t[0][-1]]-yy[t[0][0]]
print("M1: center: %f m, width: %f m, half-divergence=%f rad"%(center,width,0.5*width/radius.min()))

# print(">>>>Div M1: ",(575-75)/2000*1.605/radius.min())
# print(">>>>Half Div M1: ",0.5*(575-75)/2000*1.605/radius.min())
print(">>>>B min: ",B.min())
# plot(yy,B2,title="radius")




