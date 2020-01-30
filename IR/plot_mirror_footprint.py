
import Shadow
import numpy
from srxraylib.plot.gol import plot_image, plot_scatter, plot, plot_show

def lateral_histogram(H,coordinate='Y',show=1, filename="", title=""):

    if coordinate == 'Y':
        direction = "h"
    else:
        direction = 'v'
    # CALCULATE fwhm h
    hh = h["histogram_%s"%direction]
    bin = h["bin_%s_center"%direction]

    tt = numpy.where(hh >= max(hh) * 0.1)

    if hh[tt].size > 1:
        binSize = bin[1] - bin[0]
        fwp1m = binSize * (tt[0][-1] - tt[0][0])
        fwp1m_coordinates = (bin[tt[0][0]], bin[tt[0][-1]])
    print(fwp1m, fwp1m_coordinates)
    xt = numpy.array([1e3*bin.min(),1e3*bin.max()])
    yt = numpy.array([0.1*hh.max(),0.1*hh.max()])
    f = plot(1e3*bin, hh, xt, yt, xtitle="%s [mm]"%coordinate, title=title+"; Size at 0.1*height: %7.1f mm"%(1e3*fwp1m),show=0)
    if filename != "":
        f[0].savefig(filename)
        print("File written to disk: ",filename)
    if show: plot_show()
    return 1e3*fwp1m

if __name__ == "__main__":

    do_show = 0
    WX = []
    WY = []

    for i in range(1,6):
        filename = "/home/manuel/Oasys/mirr.%02d"%i
        print(">>>>",filename)
        beam = Shadow.Beam()
        beam.load(filename)

        h = beam.histo2(2, 1, ref=23, nbins=200)

        for key in h.keys():
            print(key)

        f = plot_image(h["histogram"], 1e3*h["bin_h_center"], 1e3*h["bin_v_center"],xtitle="Y [mm]",ytitle="X [mm]",title=filename,show=0)
        f[0].savefig("footprints/footprint_mirror%d"%i)
        if do_show:
            plot_show()

        WY.append( lateral_histogram(h,'Y',show=do_show, filename="footprints/footprint_mirror%d_histo_Y"%i, title="Mirror %d"%i) )
        WX.append( lateral_histogram(h,'X',show=do_show, filename="footprints/footprint_mirror%d_histo_X"%i, title="Mirror %d"%i) )


    for i in range(5):
        print("dimensions (width at 0.1 peak) for mirror %d: Y: %f mm X: %f mm."%(i+1,WY[i],WX[i]))
