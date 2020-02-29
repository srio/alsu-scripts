import numpy
from srxraylib.plot.gol import  plot
import matplotlib.pylab as plt

import matplotlib
matplotlib.rc('xtick', labelsize=15)
matplotlib.rc('ytick', labelsize=15)
matplotlib.rcParams.update({'font.size': 15})


filenames = ["deformationcryogenic1d","deformationcryogenic1dKh","deformationwater1_1d","deformationwater2_1d"]
rangex=[75,75,60,60]
rangey=[0.25,0.25,None,None]
# filename = "cryogenic1d.txt"

for i,filename in enumerate(filenames):

    a = numpy.loadtxt(filename+".txt")
    if rangey[i] is None:
        yrange=None
    else:
        yrange = [-rangey[i],rangey[i]]

    fig = plot(1e3*a[:,0],1e6*a[:,1],xtitle="w [mm]",ytitle="height [$\mu$m]",xrange=[-rangex[i],rangex[i]],yrange=yrange,show=0)
    fig[0].subplots_adjust(bottom=0.15)

    plt.savefig(filename+".png")
    print("File %s.png written to file"%filename)
    plt.show()