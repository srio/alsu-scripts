import Shadow

def run_shadow(monochromatic=0):
    #
    # Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
    #
    import Shadow
    import numpy

    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 1

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

    oe0.BENER = 2.0
    oe0.EPSI_X = 1e-12
    oe0.EPSI_Z = 1e-12
    oe0.FDISTR = 4
    oe0.FSOURCE_DEPTH = 4
    oe0.F_COLOR = 3
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.0005
    oe0.HDIV2 = 0.0005
    oe0.ISTAR1 = 5676561
    oe0.NCOL = 0
    oe0.NPOINT = 100000
    oe0.N_COLOR = 0
    if monochromatic:
        oe0.PH1 = 12400.000000
        oe0.PH2 = 12400.000001

    else:
        oe0.PH1 = 12395.0
        oe0.PH2 = 12405.0
    oe0.POL_DEG = 0.0
    oe0.R_ALADDIN = 133.4256380792608
    oe0.R_MAGNET = 1.3342563807926082
    oe0.SIGDIX = 0.0
    oe0.SIGDIZ = 0.0
    oe0.SIGMAX = 0.0001
    oe0.SIGMAY = 0.0
    oe0.SIGMAZ = 0.0001
    oe0.VDIV1 = 0.0005
    oe0.VDIV2 = 0.0005
    oe0.WXSOU = 0.0
    oe0.WYSOU = 0.0
    oe0.WZSOU = 0.0

    oe1.DUMMY = 1.0
    oe1.FCYL = 1
    oe1.FMIRR = 4
    oe1.FWRITE = 1
    oe1.F_DEFAULT = 0
    oe1.SIMAG = 10000000000.0
    oe1.SSOUR = 650.0
    oe1.THETA = 89.7421689922
    oe1.T_IMAGE = 0.0
    oe1.T_INCIDENCE = 89.7421689922
    oe1.T_REFLECTION = 89.7421689922
    oe1.T_SOURCE = 650.0

    oe2.DUMMY = 1.0
    oe2.FILE_REFL = b'/home/manuel/Oasys/bragg.dat'
    oe2.FWRITE = 1
    oe2.F_CENTRAL = 1
    oe2.F_CRYSTAL = 1
    oe2.PHOT_CENT = 12400.0
    oe2.R_LAMBDA = 5000.0
    oe2.T_IMAGE = 0.0
    oe2.T_INCIDENCE = 80.824080911
    oe2.T_REFLECTION = 80.824080911
    oe2.T_SOURCE = 0.0

    oe3.ALPHA = 180.0
    oe3.DUMMY = 1.0
    oe3.FILE_REFL = b'/home/manuel/Oasys/bragg.dat'
    oe3.FILE_RIP = b'//home/manuel/OASYS1.2/alsu-scripts/GEMINI/BL822_crystal_error_shadow.dat'
    oe3.FWRITE = 1
    oe3.F_CENTRAL = 1
    oe3.F_CRYSTAL = 1
    oe3.F_G_S = 2
    oe3.PHOT_CENT = 12400.0
    oe3.R_LAMBDA = 5000.0
    oe3.THICKNESS = 10.0
    oe3.T_IMAGE = 0.0
    oe3.T_INCIDENCE = 80.6539124458
    oe3.T_REFLECTION = 80.9936764183
    oe3.F_EXT = 1


    oe3.FWRITE = 0
    oe3.R_JOHANSSON = -333.33
    oe3.A_BRAGG = 0 #.17
    oe3.F_BRAGG_A = 1
    oe3.F_JOHANSSON = 1
    oe3.F_RIPPLE = 1

    # oe3.ALPHA = 180.0
    # oe3.DUMMY = 1.0
    # oe3.FILE_REFL = b'/home/manuel/Oasys/bragg.dat'
    # oe3.FWRITE = 1
    # oe3.F_CENTRAL = 1
    # oe3.F_CRYSTAL = 1
    # oe3.PHOT_CENT = 12400.0
    # oe3.R_LAMBDA = 5000.0
    # oe3.T_IMAGE = 0.0
    # oe3.T_INCIDENCE = 80.824080911
    # oe3.T_REFLECTION = 80.824080911


    oe4.DUMMY = 1.0
    oe4.FCYL = 1
    oe4.FMIRR = 4
    oe4.F_DEFAULT = 0
    oe4.F_SIDE = 1
    oe4.SIMAG = 977.0
    oe4.SSOUR = 10000000000.0
    oe4.THETA = 89.7421689922
    oe4.T_IMAGE = 977.0
    oe4.T_INCIDENCE = 89.7421689922
    oe4.T_REFLECTION = 89.7421689922
    oe4.T_SOURCE = 293.0
    oe4.FWRITE = 0

    # Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    beam.genSource(oe0)

    # import scipy.constants as codata
    #
    # energies = [12395, 12405]
    #
    #
    # nenergies = len(energies)
    #
    # for i in range(beam.nrays()):
    #     energy = energies[numpy.mod(i, nenergies)]
    #     wavelength_in_m = codata.c * codata.h / codata.e / energy
    #     beam.rays[i, 10] = 2 * numpy.pi / (wavelength_in_m * 100)



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
    import os
    os.system("rm start.* begin.dat star.* mirr.* rmirr.* optax.* effic.* end.*")
    beam = run_shadow()


    os.system("rm star.* mirr.* rmirr.* optax.* effic.* end.*")
    os.system("./shadow3/shadow3  < shadow3.inp")
    beam = Shadow.Beam()
    beam.load("star.04")





    Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, ref=23, title="Real space") #, yrange=[-60e-4,60e-4])
    # Shadow.ShadowTools.plotxy(beam, 4, 6, nbins=101, nolost=1, ref=23, title="Divergences")
    # Shadow.ShadowTools.histo1(beam,11,ref=23,nolost=1,xrange=[12395.0,12405])
