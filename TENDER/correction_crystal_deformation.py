#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow
import numpy
from srxraylib.plot.gol import set_qt, plot

set_qt()

def run_shadow(FILE_RIP="",
               SIMAG4=3.842,
               SSOUR4=26.334,
               SIMAG5=3.343,
               SSOUR5=26.834,
               ):
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
    oe5 = Shadow.OE()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.FDISTR = 3
    oe0.F_COLOR = 3
    oe0.F_PHOT = 0
    oe0.F_POLAR = 0
    oe0.HDIV1 = 0.0
    oe0.HDIV2 = 0.0
    oe0.IDO_VX = 0
    oe0.IDO_VZ = 0
    oe0.IDO_X_S = 0
    oe0.IDO_Y_S = 0
    oe0.IDO_Z_S = 0
    oe0.ISTAR1 = 123456
    oe0.NPOINT = 100000
    oe0.PH1 = 2999.5
    oe0.PH2 = 3000.5
    oe0.SIGDIX = 8e-06
    oe0.SIGDIZ = 8e-06
    oe0.SIGMAX = 1.5e-05
    oe0.SIGMAZ = 1.5e-05
    oe0.VDIV1 = 0.0
    oe0.VDIV2 = 0.0

    oe1.ALPHA = 270.0
    oe1.DUMMY = 100.0
    oe1.FWRITE = 1
    oe1.T_IMAGE = 0.0
    oe1.T_INCIDENCE = 89.5416337639
    oe1.T_REFLECTION = 89.5416337639
    oe1.T_SOURCE = 10.7324

    oe2.ALPHA = 270.0
    oe2.DUMMY = 100.0
    oe2.FHIT_C = 1
    oe2.FILE_REFL = b'C:\\Users\\Manuel\\Oasys/bragg.dat'
    oe2.FILE_RIP = FILE_RIP
    oe2.FWRITE = 1
    oe2.F_CENTRAL = 1
    oe2.F_CRYSTAL = 1
    oe2.F_G_S = 2
    oe2.F_RIPPLE = 1
    oe2.PHOT_CENT = 3000.0
    oe2.RLEN1 = 0.009937908648111332
    oe2.RLEN2 = 0.009937908648111332
    oe2.RWIDX1 = 0.009878988142292489
    oe2.RWIDX2 = 0.009878988142292489
    oe2.R_LAMBDA = 5000.0
    oe2.T_IMAGE = 0.0
    oe2.T_INCIDENCE = 75.777
    oe2.T_REFLECTION = 75.777
    oe2.T_SOURCE = 10.1196

    oe3.ALPHA = 180.0
    oe3.DUMMY = 100.0
    oe3.FILE_REFL = b'C:\\Users\\Manuel\\Oasys\\bragg.dat'
    oe3.FWRITE = 1
    oe3.F_CENTRAL = 1
    oe3.F_CRYSTAL = 1
    oe3.PHOT_CENT = 3000.0
    oe3.R_LAMBDA = 5000.0
    oe3.T_IMAGE = -0.009
    oe3.T_INCIDENCE = 75.777
    oe3.T_REFLECTION = 75.777
    oe3.T_SOURCE = 0.009

    oe4.DUMMY = 100.0
    oe4.FCYL = 1
    oe4.FMIRR = 2
    oe4.FWRITE = 1
    oe4.F_DEFAULT = 0
    oe4.SIMAG = SIMAG4
    oe4.SSOUR = SSOUR4
    oe4.THETA = 89.5416337639
    oe4.T_IMAGE = 0.0
    oe4.T_INCIDENCE = 89.5416337639
    oe4.T_REFLECTION = 89.5416337639
    oe4.T_SOURCE = 5.482

    oe5.ALPHA = 270.0
    oe5.DUMMY = 100.0
    oe5.FCYL = 1
    oe5.FMIRR = 2
    oe5.FWRITE = 1
    oe5.F_DEFAULT = 0
    oe5.SIMAG = SIMAG5
    oe5.SSOUR = SSOUR5
    oe5.THETA = 89.5416337639
    oe5.T_IMAGE = 3.343
    oe5.T_INCIDENCE = 89.5416337639
    oe5.T_REFLECTION = 89.5416337639
    oe5.T_SOURCE = 0.5

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

    #
    # run optical element 5
    #
    print("    Running optical element: %d" % (5))
    if iwrite:
        oe5.write("start.05")

    beam.traceOE(oe5, 5)

    if iwrite:
        oe5.write("end.05")
        beam.write("star.05")
    return beam

if __name__ == "__main__":

    SIMAG4 = 3.842
    SSOUR4 = 26.334
    SIMAG5 = 3.343
    SSOUR5 = 26.834

    FILE_DIR = 'C:/Users/Manuel/Oasys/TENDER DCM Performance/'
    FILE_RIP = 'disp1000_shadow.dat'

    FILE_RIP_LIST = []
    FILE_RIP_LIST.append('disp100000_shadow.dat')
    FILE_RIP_LIST.append('disp10000_shadow.dat')
    FILE_RIP_LIST.append('disp3000_shadow.dat')
    FILE_RIP_LIST.append('disp2500_shadow.dat')
    FILE_RIP_LIST.append('disp2000_shadow.dat')
    FILE_RIP_LIST.append('disp1500_shadow.dat')
    FILE_RIP_LIST.append('disp1250_shadow.dat')
    FILE_RIP_LIST.append('disp1000_shadow.dat')

    nruns = 20


    # defocus = numpy.linspace(-0.25, 0.50, nruns)
    defocus = numpy.linspace(-0.75, 0.50, nruns)
    # defocus = numpy.linspace(-1, 1, nruns)

    f = open('correction_crystal_deformation.txt','w')

    for FILE_RIP in FILE_RIP_LIST:
        out = numpy.zeros((nruns, 4))
        for i in range(nruns):
            beam = run_shadow(FILE_RIP=(FILE_DIR+FILE_RIP).encode(),
                            SIMAG4=  3.842  + defocus[i],
                            SSOUR4 = 26.334,
                            SIMAG5 = 3.343  + defocus[i],
                            SSOUR5 = 26.834,
                            )

            tkt = beam.histo2(3, 1, nbins=100, nolost=1, ref=23, )

            # for key in tkt.keys():
            #     print(key)

            print("FWHM: %f um (H) x %f um (V)" % (1e6 * tkt["fwhm_h"], 1e6 * tkt["fwhm_v"]))

            out[i, 0] = defocus[i]
            out[i, 1] = defocus[i]
            out[i, 2] = 1e6 * tkt["fwhm_h"]
            out[i, 3] = 1e6 * tkt["fwhm_v"]

            # Shadow.ShadowTools.plotxy(beam, 3, 1, nbins=101, nolost=1, title="Real space")

        # plot(out[:, 0], out[:, 2], out[:, 0], out[:, 3])

        i_opt = numpy.argmin(out[:, 2])
        j_opt = numpy.argmin(out[:, 3])

        f.write("File: %s \n" % (FILE_RIP))
        f.write("Optimum at defocus: %f m (H) x %f m (V)\n" % (out[i_opt, 0], out[j_opt, 0]))
        f.write("Optimum FWHM: %f m (H) x %f m (V)\n\n" % (out[i_opt, 2], out[j_opt, 3]))


    f.close()
    print("File correction_crystal_deformation.txt written to disk.")