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

y -= y[y.size//2]
y *= 1e-3

plot(y, B, y, B2, legend=["original","smoothed"],xtitle="y / m",ytitle="B / T")

f = open("BM_multi.b", "w")
for i in range(y.size):
    f.write("%f  %f\n" % (y[i], B2[i]))
f.close()
print("File written to disk: BM_multi.b")


