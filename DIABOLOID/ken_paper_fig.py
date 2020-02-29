


import Shadow
from srxraylib.plot.gol import plot_scatter
import matplotlib.pylab as plt
import numpy
from matplotlib import rc

# # rc('font', **{'family':'serif','serif':['Palatino']})
# rc('text', usetex=True)



b = Shadow.Beam()
b.load("/home/manuel/Oasys/star_ken.01")

x = b.getshonecol(1,nolost=1)
y = b.getshonecol(3,nolost=1)

f = plot_scatter(1e6*x,1e6*y,xrange=[-500,500],yrange=[-500,500],nbins=51,xtitle='X [$\mu$ m]',ytitle='y [$\mu$ m]',show=0)


# change in gol.py
# if plot_histograms:
#     left, width = 0.12, 0.65
#     bottom, height = 0.1, 0.65


f[1].set_ylabel("y [$\mu$ m]",fontsize=15)
f[1].set_xlabel("x [$\mu$ m]",fontsize=15)
f[1].tick_params(labelsize=15)
f[2].tick_params(labelsize=0)
f[3].tick_params(labelsize=0)

f[0].savefig("Diaboloid_raytrace1.png",dpi=600)
# f[0].savefig("raytracing.pdf",dpi=600)


plt.show()
