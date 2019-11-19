

import Shadow
import numpy
from srxraylib.plot.gol import plot, plot_image

beam = Shadow.Beam()

beam.load("/home/manuel/Oasys/star.02")

# Shadow.ShadowTools.plotxy(beam,1,3,xrange=[-100e-6,100e-6],yrange=[-4000e-6,4000e-6],nbins=300)


nbins=1000
npositions = 300
max_x = 5
positions = numpy.linspace(-max_x,max_x,npositions)

out = numpy.zeros((nbins,npositions))



for i in range(npositions):
    beami = beam.duplicate()
    beami.retrace(positions[i])
    tkt = beami.histo1(1,xrange=[-0.28,0.28],nbins=nbins,nolost=1,ref=23)
    out[:,i] = tkt["histogram"]


for t in tkt.keys():
    print(t)

# beam = Shadow.Beam()
# beam.load("/home/manuel/Oasys/star.02")
# tkt = beam.histo1(1,xrange=[-80000e-6,5000e-6],nbins=nbins,nolost=1,ref=23)
# Shadow.ShadowTools.histo1(beam,1,nbins=nbins,nolost=1,ref=23)
# plot(tkt["bin_center"],tkt["histogram"])



plot_image(out,tkt["bin_center"],positions,aspect='auto',xtitle="Horizontal X [m]",ytitle="Propagation distance [m]",show=False,figsize=[10,2])
plot_image(numpy.log10(out),tkt["bin_center"],positions,aspect='auto',xtitle="Horizontal X [m]",ytitle="Propagation distance [m]",title="log",show=True)


