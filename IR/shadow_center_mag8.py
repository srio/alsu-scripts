import Shadow
import numpy

def get_R(p_foc,q_foc,incidence):
    mm = (1.0 / p_foc + 1.0 / q_foc)
    return 2 / (numpy.cos(incidence * numpy.pi / 180)) / mm

def get_q(R,p_foc,incidence):
    qq = 2 / (R * numpy.cos(incidence * numpy.pi / 180))  - 1.0 / p_foc
    return 1.0 / qq

def run_shadow(angle=180.0):
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

    # oe0.BENER = 2.0
    # oe0.EPSI_X = 7e-11
    # oe0.EPSI_Z = 7e-11
    # oe0.FDISTR = 6
    # oe0.FSOURCE_DEPTH = 4
    # oe0.F_COLOR = 3
    # oe0.F_PHOT = 0
    # oe0.HDIV1 = 0.033
    # oe0.HDIV2 = 0.033
    # oe0.ISTAR1 = 5676561
    # oe0.NCOL = 0
    # oe0.NPOINT = 300000
    # oe0.N_COLOR = 0
    # oe0.PH1 = 0.4
    # oe0.PH2 = 0.401
    # oe0.POL_DEG = 0.0
    # oe0.R_ALADDIN = -7.615618611829955
    # oe0.R_MAGNET = -7.615618611829955
    # oe0.SIGDIX = 0.0
    # oe0.SIGDIZ = 0.0
    # oe0.SIGMAX = 7e-06
    # oe0.SIGMAY = 0.0
    # oe0.SIGMAZ = 1e-05
    # oe0.VDIV1 = 0.05
    # oe0.VDIV2 = 0.05
    # oe0.WXSOU = 0.0
    # oe0.WYSOU = 0.0
    # oe0.WZSOU = 0.0

    oe1.ALPHA = 90.0
    oe1.DUMMY = 100.0
    oe1.FMIRR = 3
    oe1.FWRITE = 1
    oe1.T_IMAGE = 2.935
    oe1.T_INCIDENCE = 45.0
    oe1.T_REFLECTION = 45.0
    oe1.T_SOURCE = 2.935

    oe2.ALPHA = 90.0
    oe2.DUMMY = 100.0
    oe2.FWRITE = 1
    oe2.F_REFRAC = 2
    oe2.T_IMAGE = 0.0
    oe2.T_INCIDENCE = 0.0
    oe2.T_REFLECTION = 180.0
    oe2.T_SOURCE = 0.0

    # oe3.DUMMY = 100.0
    # oe3.FWRITE = 1
    # oe3.F_REFRAC = 2
    # oe3.T_IMAGE = 1.0
    # oe3.T_INCIDENCE = 0.0
    # oe3.T_REFLECTION = angle
    # oe3.T_SOURCE = 0.0


    #
    # last arm
    #
    oe3.ALPHA = 90.0
    oe3.DUMMY = 100.0
    oe3.FWRITE = 1
    oe3.F_REFRAC = 2
    oe3.T_IMAGE = 0.0
    oe3.T_INCIDENCE = 0.0
    oe3.T_REFLECTION = angle
    oe3.T_SOURCE = 0.0

    oe4.ALPHA = 270.0
    oe4.DUMMY = 100.0
    oe4.FWRITE = 1
    oe4.F_REFRAC = 2
    oe4.T_IMAGE = 0.0
    oe4.T_INCIDENCE = 0.0
    oe4.T_REFLECTION = 180.0
    oe4.T_SOURCE = 0.0



    # Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    beam.load("/home/manuel/Oasys/begin_ir_0p1_100000.dat")

    y = beam.rays[:, 1]
    select_mode = 2  # 0=all, 1=Mag7, 2=Mag8, 3=both RB only, 4=RB1, 5=RB2

    if select_mode == 0:
        pass
    elif select_mode == 1:
        ibad = numpy.where(y < -0.3)
        beam.rays[ibad, 6:9] = 0.0
        beam.rays[ibad, 15:18] = 0.0

        ibad = numpy.where(y > 0.3)
        beam.rays[ibad, 6:9] = 0.0
        beam.rays[ibad, 15:18] = 0.0

    elif select_mode == 2:
        ibad = numpy.where(y < 0.66)
        beam.rays[ibad, 6:9] = 0.0
        beam.rays[ibad, 15:18] = 0.0
    elif select_mode == 3:
        ibad = numpy.where(numpy.abs(y) < 0.3)
        beam.rays[ibad, 6:9] = 0.0
        beam.rays[ibad, 15:18] = 0.0

        ibad = numpy.where(y > 0.66)
        beam.rays[ibad, 6:9] = 0.0
        beam.rays[ibad, 15:18] = 0.0
    elif select_mode == 4:
        ibad = numpy.where(y > -0.3)
        beam.rays[ibad, 6:9] = 0.0
        beam.rays[ibad, 15:18] = 0.0

    elif select_mode == 5:
        ibad = numpy.where(y < 0.3)
        beam.rays[ibad, 6:9] = 0.0
        beam.rays[ibad, 15:18] = 0.0

        ibad = numpy.where(y > 0.66)
        beam.rays[ibad, 6:9] = 0.0
        beam.rays[ibad, 15:18] = 0.0


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
    beam.traceOE(oe4, 4)

    if iwrite:
        oe3.write("end.04")
        beam.write("star.04")



    # Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space")
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")



    return beam


def find_q_for_Mag8():

    angles = numpy.loadtxt("/home/manuel/Oasys/angle.01")
    print(">>>>",angles.shape)
    # plot(angles[:,0],angles[:,1])
    print("Mean incident angle: ", angles[:, 1].mean())
    print("Mean ouput angle: ", angles[:, 2].mean())

    p0 = 2.935
    q0 = p0
    p1 = p0 + 0.15
    q1 = q0 - 0.15
    R = get_R(p0,q0,45.0)
    print("R toroid:",)
    print("q toroid:", get_q(R, p0, 45.0))
    print("mean angle ",angles[:, 1].mean())

    q1 = get_q(R, p1, angles[:, 1].mean())
    print("q focus Mag8:", q1 )
    return q1


if __name__ == "__main__":

    import Shadow
    from srxraylib.plot.gol import plot

    nruns = 30
    angle = numpy.linspace(175,185,nruns)
    xc = numpy.zeros_like(angle)
    zc = numpy.zeros_like(angle)

    for i in range(nruns):
        beam = run_shadow(angle=angle[i])
        beam.retrace(10)
        xc[i] = beam.rays[:, 0].mean()
        zc[i] = beam.rays[:, 3].mean()


    plot(angle, xc, show=0, xtitle="X")
    plot(angle, zc, show=1, xtitle="Z")

    # print("Distance: ",find_q_for_Mag8() - 2.935)


