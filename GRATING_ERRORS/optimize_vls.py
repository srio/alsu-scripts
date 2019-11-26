
import numpy

def run_shadow(b2=0.0):
    #
    # Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
    #
    import Shadow
    import numpy

    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 0

    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = Shadow.Beam()
    oe0 = Shadow.Source()
    oe1 = Shadow.OE()
    oe2 = Shadow.OE()
    oe3 = Shadow.OE()
    oe4 = Shadow.OE()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.FDISTR = 3
    oe0.F_COHER = 1
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.0
    oe0.HDIV2 = 0.0
    oe0.IDO_VX = 0
    oe0.IDO_VZ = 0
    oe0.IDO_X_S = 0
    oe0.IDO_Y_S = 0
    oe0.IDO_Z_S = 0
    oe0.ISTAR1 = -2085389133
    oe0.NPOINT = 50000
    oe0.PH1 = 500.0
    oe0.SIGDIX = 1.6e-05
    oe0.SIGDIZ = 1.6e-05
    oe0.SIGMAX = 0.0176
    oe0.SIGMAZ = 0.0176
    oe0.VDIV1 = 0.0
    oe0.VDIV2 = 0.0

    oe1.DUMMY = 0.1
    oe1.FCYL = 1
    oe1.FHIT_C = 1
    oe1.FMIRR = 2
    oe1.FWRITE = 1
    oe1.F_DEFAULT = 0
    oe1.RLEN1 = 500.0
    oe1.RLEN2 = 500.0
    oe1.RWIDX1 = 100.0
    oe1.RWIDX2 = 100.0
    oe1.SIMAG = 14100.0
    oe1.SSOUR = 110400.0
    oe1.THETA = 88.951
    oe1.T_IMAGE = 5000.0
    oe1.T_INCIDENCE = 88.951
    oe1.T_REFLECTION = 88.951
    oe1.T_SOURCE = 110400.0

    oe2.DUMMY = 0.1
    oe2.FWRITE = 1
    oe2.T_IMAGE = 600.0
    oe2.T_INCIDENCE = 88.775
    oe2.T_REFLECTION = 88.775
    oe2.T_SOURCE = 250.0

    oe3.DUMMY = 0.1
    oe3.FWRITE = 3
    oe3.F_REFRAC = 2
    oe3.F_SCREEN = 1
    oe3.N_SCREEN = 1
    oe3.T_IMAGE = 0.0
    oe3.T_INCIDENCE = 0.0
    oe3.T_REFLECTION = 180.0
    oe3.T_SOURCE = 0.0

    oe4.DUMMY = 0.1
    oe4.FHIT_C = 1
    oe4.FWRITE = 1
    oe4.F_GRATING = 1
    oe4.F_RULING = 5
    oe4.F_RUL_ABS = 1
    oe4.RLEN1 = 120.0
    oe4.RLEN2 = 120.0
    oe4.RULING = 300.0
    oe4.RUL_A1 = 0.0257
    oe4.RUL_A2 = b2
    oe4.RWIDX1 = 25.0
    oe4.RWIDX2 = 25.0
    oe4.T_IMAGE = 19600.0
    oe4.T_INCIDENCE = 89.273191
    oe4.T_REFLECTION = 87.673376
    oe4.T_SOURCE = 175.0

    # Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")

    #
    # run optical element 1
    #
    print("    Running optical element: %d" % (1))
    if iwrite:
        oe1.write("start.01")

    beam.traceOE(oe1, 1)

    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")

    #
    # run optical element 2
    #
    print("    Running optical element: %d" % (2))
    if iwrite:
        oe2.write("start.02")

    beam.traceOE(oe2, 2)

    if iwrite:
        oe2.write("end.02")
        beam.write("star.02")

    #
    # run optical element 3
    #
    print("    Running optical element: %d" % (3))
    if iwrite:
        oe3.write("start.03")

    beam.traceOE(oe3, 3)

    if iwrite:
        oe3.write("end.03")
        beam.write("star.03")

    #
    # run optical element 4
    #
    print("    Running optical element: %d" % (4))
    if iwrite:
        oe4.write("start.04")

    beam.traceOE(oe4, 4)

    if iwrite:
        oe4.write("end.04")
        beam.write("star.04")

    # Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space")
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")
    return beam

if __name__ == "__main__":

    from srxraylib.plot.gol import plot
    import Shadow


    B2 = numpy.linspace(-1e-5, 1e-5, 100)
    OUT = numpy.zeros_like(B2)



    for i,b2 in enumerate(B2):

        beam = run_shadow(b2)
        # Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space b2: %f"%b2)
        col3 = beam.getshonecol(3,nolost=1)
        OUT[i] = col3.std()


        #
        # print("\n\nset energy: %f eV"%energy1)
        # print("alpha: ",alpha1)
        # print("beta: ",beta1)
        #
        # dict = respower(beam, 11, 1, hlimit=0.1, nolost=True)
        # dict5 = respower(beam, 11, 1, hlimit=0.5, nolost=True)
        #
        #
        # # store results
        # Energy[i] = energy1
        # Resolution[i] = dict["resolvingPower"]
        # Resolution5[i] = dict5["resolvingPower"]
        # Alpha[i] = alpha1
        # Beta[i] = beta1
        # C[i] = numpy.cos(beta1*numpy.pi/180) / numpy.cos(alpha1*numpy.pi/180)
        # M[i] = 7.573 / 25.201 / C[i]
        #
        #
        # # run now the monochromatic case to get the monochromatic size
        # # beam1 = run_shadow(energy1, 0.0,  alpha1, beta1, sigmas)
        # # tkt = beam1.histo1(1,nolost=True, ref=23, nbins=100)
        # tkt = dict["histo_dict"]
        # # plot(tkt["bin_path"] * 1e6, tkt["histogram_path"], title="E=%s  S=%s"%(energy1,tkt["fwhm"] * 1e6))
        # # ImageSize[i] = tkt["fwhm_subpixel"]*1e6
        # ImageSize[i] = tkt["fwhm"]*1e6
        # SourceSize[i] = 2.35 * sigmas[1] * 1e6
        #
        # # get now the footprint on grating
        # b = Shadow.Beam()
        # b.load("mirr.03")
        # tkt_m = b.histo1(2, nolost=True)
        # # Footprint[i] = tkt_m["fwhm_subpixel"]
        # Footprint[i] = tkt_m["fwhm"]


    plot(B2,OUT)
    # #
    # # dump to file
    # #
    #
    # f = open(dumpfile,'w')
    # f.write("\n#S 1\n")
    # f.write("#N 10\n")
    # f.write("#L photon energy [eV]  alpha [deg]  beta [deg]  size source V [um]  size [um]  source*M  footprint on grating [m]  C  resolving power 0.1  resolving power 0.5\n")
    # for i in range(Energy.size):
    #     f.write("%f   %f   %f   %f   %f   %f   %f   %f   %f   %f\n"%(
    #         Energy[i],
    #         Alpha[i],
    #         Beta[i],
    #         SourceSize[i],
    #         ImageSize[i],
    #         SourceSize[i]*M[i],
    #         Footprint[i],
    #         C[i],
    #         Resolution[i],
    #         Resolution5[i]))
    # f.close()
    # print("File written to disk: %s"%dumpfile)