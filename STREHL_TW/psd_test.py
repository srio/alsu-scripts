import numpy
from srxraylib.metrology.dabam import dabam
from srxraylib.plot.gol import plot

def fitpowerlaw(f, psdHeights):
    # fit PSD to a power law
    x = numpy.log10(f)
    y_h = numpy.log10(psdHeights)
    # y_s = numpy.log10(psdSlopes)
    # select the fitting area (80% of the full interval, centered)

    x_left  = x.min() #(x.min() + 0.1 * (x.max() - x.min()))
    x_right = x.max() #(x.max() - 0.1 * (x.max() - x.min()))

    # redefine  left limit for the fit to the frequency value corresponding to the correlation length
    # acf_h = autocovariance_1D(self.sy,self.zprof)
    # f_h = numpy.log10( 1.0 / acf_h[2] )
    # x_left = f_h

    c1 = (x < x_right)
    c2 = (x > x_left)
    igood = numpy.where(c1 & c2)
    igood = numpy.array(igood)
    igood.shape = -1

    coeffs_h = numpy.polyfit(x[igood], y_h[igood], 1)
    # coeffs_s = numpy.polyfit(x[igood], y_s[igood], 1)

    powerlaw = {"hgt_pendent": coeffs_h[0], "hgt_shift": coeffs_h[1], \
                     "index_from": igood[0], "index_to": igood[-1]}

    return powerlaw

def psd(xx, yy, onlyrange = None):
    """
     psd: Calculates the PSD (power spectral density) from a profile

      INPUTS:
           x - 1D array of (equally-spaced) lengths.
           y - 1D array of heights.
      OUTPUTS:
           f - 1D array of spatial frequencies, in units of 1/[x].
           s - 1D array of PSD values, in units of [y]^3.
      KEYWORD PARAMETERS:
           onlyrange - 2-element array specifying the min and max spatial
               frequencies to be considered. Default is from
               1/(length) to 1/(2*interval) (i.e., the Nyquist
               frequency), where length is the length of the scan,
               and interval is the spacing between points.

      PROCEDURE
            Use FFT

    """
    n_pts = xx.size
    if (n_pts <= 1):
        print ("psd: Error, must have at least 2 points.")
        return 0

    window=yy*0+1.
    length=xx.max()-xx.min()  # total scan length.
    delta = xx[1] - xx[0]

    # psd from windt code
    # s=length*numpy.absolute(numpy.fft.ifft(yy*window)**2)
    # s=s[0:(n_pts/2+1*numpy.mod(n_pts,2))]  # take an odd number of points.

    #xianbo + luca:
    s0 = numpy.absolute(numpy.fft.fft(yy*window))
    s =  2 * delta * s0[0:int(len(s0)/2)]**2/s0.size # uniformed with IGOR, FFT is not symmetric around 0
    s[0] /= 2
    s[-1] /= 2


    n_ps=s.size                       # number of psd points.
    interval=length/(n_pts-1)         # sampling interval.
    f_min=1./length                   # minimum spatial frequency.
    f_max=1./(2.*interval)            # maximum (Nyquist) spatial frequency.
    # spatial frequencies.
    f=numpy.arange(float(n_ps))/(n_ps-1)*(f_max-f_min)+f_min

    if onlyrange != None:
        roi =  (f <= onlyrange[1]) * (f >= onlyrange[0])
        if roi.sum() > 0:
            roi = roi.nonzero()
            f = f[roi]
            s = s[roi]

    return s,f


if __name__ == "__main__":
    from srxraylib.plot.gol import plot
    import matplotlib.pylab as plt

    # dm = dabam.initialize_from_external_data("fractalprofile.dat",
    #                                          column_index_abscissas=0,
    #                                          column_index_ordinates=1,
    #                                          skiprows=0,
    #                                          useHeightsOrSlopes=0,
    #                                          to_SI_abscissas=1.0,
    #                                          to_SI_ordinates=1.0,
    #                                          detrending_flag=-1)

    # dm.plot("psd_h")

    a = numpy.loadtxt("fractalprofile.dat")
    x = a[:, 0].copy()
    y = a[:, 1].copy()
    print(x.shape, y.shape)
    # plot(x,y)

    psdHeights, f = psd(x, y)
    print(f.shape, psdHeights.shape)

    pw = fitpowerlaw(f, psdHeights)
    for key in pw.keys():
        print("%s = "%key,pw[key])

    ff = plt.figure()
    plt.loglog(f, psdHeights)
    y = (f ** pw["hgt_pendent"]) * (10 ** pw["hgt_shift"])
    i0 = pw["index_from"]
    i1 = pw["index_to"]
    plt.loglog(f, y)
    plt.loglog(f[i0:i1], y[i0:i1])
    beta = -pw["hgt_pendent"]
    plt.title("PSD of heights profile (beta=%.2f,Df=%.2f)" % (beta, (5 - beta) / 2))
    plt.xlabel("f [m^-1]")
    plt.ylabel("PSD [m^3]")
    plt.plot()
    plt.show()


    # print(s.shape, f.shape)
    #
    # plot(f, psdHeights, xlog=1, ylog=1, yrange=[1e-25, 1e-19])
    # filename = "deformation%02d.dat" % entry_number
    # f = open(filename, "w")
    # for i in range(x.size):
    #     f.write("%g %g\n" % (x[i], y[i]))
    # f.close()
    # print("File written to disk: %s" % filename)

