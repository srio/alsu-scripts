import numpy
from srxraylib.plot.gol import  plot
import matplotlib.pylab as plt
from oasys.util.oasys_util import get_fwhm

import matplotlib
matplotlib.rc('xtick',         labelsize=20)
matplotlib.rc('ytick',         labelsize=20)
matplotlib.rcParams.update({'font.size': 20})


is_fit = 1

if is_fit == 0:
    filenames = ["intensityundeformed", "intensitycryogenic", "intensitycryogenicKh",
                 "intensitywater1", "intensitywater2"]
    filepng = "intensityundeformed.png"
elif is_fit == 1:
    filenames = ["intensityundeformed","intensitycryogeniccorrected","intensitycryogenicKhcorrected",
                 "intensitywater1corrected","intensitywater2corrected"]
    filepng = "intensitycorrected.png"
elif is_fit == 2:
    filenames = ["intensityundeformed", "intensitycryogeniccorrectedfit", "intensitycryogenicKhcorrectedfit",
                 "intensitywater1correctedfit", "intensitywater2correctedfit"]
    filepng = "intensitycorrectedfit.png"
elif is_fit == 3:
    filenames = ["intensityundeformed", "intensitycryogeniccorrectedfitcrop", "intensitycryogenicKhcorrectedfitcrop",
                 "intensitywater1correctedfitcrop", "intensitywater2correctedfitcrop"]
    filepng = "intensitycorrectedfitcrop.png"

# filename = "cryogenic1d.txt"

FWHM = []
STREHL = []
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

    print(">>>>>>>",get_fwhm(a[:,1],a[:,0],))
    X.append(a[:,0])
    Y.append(a[:,1])
    tmp = get_fwhm(a[:, 1], a[:, 0], )
    FWHM.append(tmp[0])
    STREHL.append(a[:,1].max())
    # plt.savefig(filename+".png")
    # print("File %s.png written to file"%filename)
    # plt.show()

TITLES = ["undeformed","cryo H","cryo V","water H","water V"]
for i in range(1,len(filenames)):
    print("%15s    FWHM: %f    STREHL: %f"%(TITLES[i],FWHM[i],STREHL[i]/STREHL[0]))
    LEGEND.append("%s FWHM:%2.1f SR:%2.1f "%(TITLES[i],FWHM[i],STREHL[i]/STREHL[0]))


fig = plt.figure(figsize=(16,8))

plt.plot(X[1], Y[1]+150,label=LEGEND[1-1])
plt.plot(X[2], Y[2]+100,label=LEGEND[2-1])
plt.plot(X[3], Y[3]+50,label=LEGEND[3-1])
plt.plot(X[4], Y[4]+0,label=LEGEND[4-1])
# plt.title(title)
plt.xlim(-10,10)
plt.ylim(0,370)
plt.xlabel("X [$\mu$m]")
plt.ylabel("intensity [a.u.]")

# plt.subplots_adjust(bottom=0.15)
    #
    # xtitle="X [$\mu$m]", ytitle="intensity [a.u.]",legend=LEGEND,
    #
    # )
ax = plt.subplot(111)
ax.legend(bbox_to_anchor=[.7,.55])


plt.savefig(filepng)
plt.show()