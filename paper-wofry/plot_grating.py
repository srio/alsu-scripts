import numpy
from srxraylib.plot.gol import  plot
import matplotlib.pylab as plt
from oasys.util.oasys_util import get_fwhm

import matplotlib
matplotlib.rc('xtick',         labelsize=20)
matplotlib.rc('ytick',         labelsize=20)
matplotlib.rcParams.update({'font.size': 20})

is_fit = True

filename = "/home/manuel/Oasys/VLS_FLEXON_to_plot.txt"

a = numpy.loadtxt(filename)
fig = plt.figure(figsize=(16,8))
plt.plot(1e3*a[:,0], 1e9*a[:,1],marker="o" )
plt.xlabel("w [mm]")
plt.ylabel("height [nm]")
ax = plt.gca()
ax.ticklabel_format(useOffset=False)
plt.savefig("grating.png")
plt.show()


filename = "/home/manuel/Oasys/intensitygrating.txt"

a = numpy.loadtxt(filename)
fig = plt.figure(figsize=(16,8))
plt.plot(a[:,0], a[:,1] )
plt.xlabel("X [$\mu$m]")
plt.ylabel("intensity [a.u.]")
ax = plt.gca()
ax.ticklabel_format(useOffset=False)
plt.savefig("intensitygrating.png")
plt.show()
