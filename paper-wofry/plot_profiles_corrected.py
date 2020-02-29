import numpy
from srxraylib.plot.gol import  plot
import matplotlib.pylab as plt
from oasys.util.oasys_util import get_fwhm

import matplotlib
matplotlib.rc('xtick',         labelsize=20)
matplotlib.rc('ytick',         labelsize=20)
matplotlib.rcParams.update({'font.size': 20})

is_fit = True

if is_fit:
    filenames = ["correctioncryogenicfit","correctioncryogenicKhfit","correctionwater1fit","correctionwater2fit"]
else:
    filenames = ["correctioncryogenic","correctioncryogenicKh","correctionwater1","correctionwater2"]

# filename = "cryogenic1d.txt"


X = []
Y = []
LEGEND = []

for i,filename in enumerate(filenames):

    a = numpy.loadtxt(filename+".txt")
    # if rangey[i] is None:
    #     yrange=None
    # else:
    #     yrange = [-rangey[i],rangey[i]]

    # fig = plot(a[:,0],a[:,1],xtitle="X [$\mu$m]",ytitle="intensity [a.u.]",show=0) #,xrange=[-rangex[i],rangex[i]],yrange=yrange)
    # fig[0].subplots_adjust(bottom=0.15)


    X.append(a[:,0])
    Y.append(a[:,1])

TITLES = ["cryo H","cryo V","water H","water V"]
for i in range(len(filenames)):
    # print("%15s    FWHM: %f    STREHL: %f"%(TITLES[i],FWHM[i],STREHL[i]/STREHL[0]))
    LEGEND.append("%s "%(TITLES[i]))


fig = plt.figure(figsize=(16,8))

plt.plot(1e3*X[1-1], 1e6*Y[1-1]+0  ,label=LEGEND[1-1])
plt.plot(1e3*X[2-1], 1e6*Y[2-1]+0  ,label=LEGEND[2-1])
plt.plot(1e3*X[3-1], 1e6*Y[3-1]+0  ,label=LEGEND[3-1])
plt.plot(1e3*X[4-1], 1e6*Y[4-1]+0  ,label=LEGEND[4-1])
# plt.title(title)
plt.xlabel("w [mm]")
plt.ylabel("height [$\mu$m]")

# plt.subplots_adjust(bottom=0.15)
    #
    # xtitle="X [$\mu$m]", ytitle="intensity [a.u.]",legend=LEGEND,
    #
    # )
ax = plt.subplot(111)
ax.legend(bbox_to_anchor=[.49,.95])

if is_fit:
    plt.savefig("correctedprofilesfit.png")
else:
    plt.savefig("correctedprofiles.png")

plt.show()