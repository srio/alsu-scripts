#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow
import numpy

def get_R(p_foc,q_foc,incidence):
    mm = (1.0 / p_foc + 1.0 / q_foc)
    return 2 / (numpy.cos(incidence * numpy.pi / 180)) / mm


def run_shadow(incidence=67.7,radius=1.0,p0=1.0,q0=1.0):
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

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.BENER = 2.734
    oe0.EPSI_X = 6.310000000000001e-10
    oe0.EPSI_Z = 3.24e-10
    oe0.FDISTR = 6
    oe0.FSOURCE_DEPTH = 4
    oe0.F_COLOR = 3
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.04
    oe0.HDIV2 = 0.04
    oe0.ISTAR1 = 5676561
    oe0.NCOL = 0
    oe0.NPOINT = 150000
    oe0.N_COLOR = 0
    oe0.PH1 = 0.0124
    oe0.PH2 = 0.01241
    oe0.POL_DEG = 0.0
    oe0.R_ALADDIN = 5.28
    oe0.R_MAGNET = 5.28
    oe0.SIGDIX = 0.0
    oe0.SIGDIZ = 0.0
    oe0.SIGMAX = 6.31e-05
    oe0.SIGMAY = 0.0
    oe0.SIGMAZ = 3.24e-05
    oe0.VDIV1 = 0.02
    oe0.VDIV2 = 0.02
    oe0.WXSOU = 0.0
    oe0.WYSOU = 0.0
    oe0.WZSOU = 0.0

    oe1.DUMMY = 100.0
    oe1.FCYL = 1
    oe1.FHIT_C = 1
    oe1.FMIRR = 1
    oe1.FWRITE = 1
    oe1.F_DEFAULT = 0
    oe1.RLEN1 = 0.1655
    oe1.RLEN2 = 0.1655
    oe1.RWIDX1 = 0.02352
    oe1.RWIDX2 = 0.02352
    oe1.SIMAG = 6.0
    oe1.SSOUR = 6.0
    oe1.THETA = 45.0
    oe1.T_IMAGE = 0.0
    oe1.T_INCIDENCE = 45.0
    oe1.T_REFLECTION = 45.0
    oe1.T_SOURCE = 6.0

    oe2.ALPHA = 90.0
    oe2.DUMMY = 100.0
    oe2.FCYL = 1
    oe2.FHIT_C = 1
    oe2.FMIRR = 1
    oe2.FWRITE = 1
    oe2.F_EXT = 1
    oe2.RLEN1 = 0.313
    oe2.RLEN2 = 0.313
    oe2.RMIRR = radius
    oe2.RWIDX1 = 0.081
    oe2.RWIDX2 = 0.081
    oe2.T_IMAGE = q0
    oe2.T_INCIDENCE = incidence
    oe2.T_REFLECTION = incidence
    oe2.T_SOURCE = p0 - 6.0

    oe3.ALPHA = 270.0
    oe3.DUMMY = 100.0
    oe3.FWRITE = 3
    oe3.F_REFRAC = 2
    oe3.T_IMAGE = 0.0
    oe3.T_INCIDENCE = 0.0
    oe3.T_REFLECTION = 180.0
    oe3.T_SOURCE = 0.0

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

    # Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space")
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")
    return beam,oe2


if __name__ == "__main__":

    from srxraylib.plot.gol import plot, set_qt
    from srxraylib.util.h5_simple_writer import H5SimpleWriter

    set_qt()
    uvalue = numpy.linspace(0.2, 3, 101)
    Radius = numpy.zeros_like(uvalue)
    Fwhm = numpy.zeros_like(uvalue)
    Std = numpy.zeros_like(uvalue)

    dump_file = False

    if dump_file:
        h = H5SimpleWriter.initialize_file("IR_BM_shadow3_scanMagnification_moreno.h5")

    grazing = 80.27
    incidence = 90.0 - grazing

    for i in range(uvalue.size):
        p0 = 12.0 * uvalue[i] / (1 + uvalue[i])
        q0 = 12.0 - p0

        print(">>>>>>>>>>>>>>>>>>>>>>>>>",uvalue[i],p0,q0,p0+q0)
        radius = get_R(p0,q0,incidence)
        beam,oe2 = run_shadow(incidence=incidence,radius=radius,p0=p0,q0=q0)
        # Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space")

        tkt = beam.histo1(1,nbins=201,nolost=1) # xrange=[-4000e-6,4000e-6],
        print(incidence,oe2.RMIRR,1e6*tkt["fwhm"])
        Radius[i] = oe2.RMIRR
        Fwhm[i] = 1e6*tkt["fwhm"]
        Std[i] = 1e6 * beam.get_standard_deviation(1,nolost=1)

        if dump_file:
            h.create_entry("iteration %f" % uvalue[i], nx_default="histogram")
            h.add_dataset(1e6*tkt["bin_path"],tkt["histogram_path"], dataset_name="histogram",entry_name="iteration %f" % uvalue[i],
                          title_x="X / um",title_y="intensity / a.u.")



    plot(uvalue,Fwhm,xtitle="uvalue",ytitle="FWHM/um",show=False)
    plot(uvalue, Fwhm,uvalue, Std,xtitle="uvalue",ytitle="FWHM/um",legend=["fwhm","std"],show=False)
    plot(uvalue, Fwhm*uvalue, uvalue, Std*uvalue, xtitle="uvalue", ytitle="FWHM/um times uvalue", legend=["fwhm", "std"])

    if dump_file:
        h.create_entry("scan results", nx_default="Fwhm")
        h.add_dataset(uvalue, Fwhm, dataset_name="Fwhm", entry_name="scan results",
                      title_x="uvalue", title_y="fwhm / um")
        h.add_dataset(uvalue, Std, dataset_name="Std", entry_name="scan results",
                      title_x="uvalue", title_y="std / um")
        h.add_dataset(Radius, Fwhm, dataset_name="Radius", entry_name="scan results",
                      title_x="radius / m", title_y="fwhm / um")

    for key in tkt.keys():
        print(key)
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")
    imin = Std.argmin()
    optimized_uvalue = uvalue[imin]     #
    p0 = 12.0 * optimized_uvalue / (1 + optimized_uvalue)
    q0 = 12.0 - p0
    radius = get_R(p0, q0, incidence)

    beam,oe1 = run_shadow(incidence=incidence,radius=radius,p0=p0,q0=q0)


    tkt = Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=201, nolost=1, ref=23,
                              # xrange=[-5e-3*5,5e-3*5],
                              title="theta: %f, R: %f, uvalue: %f"%(incidence,radius,optimized_uvalue))

    print("best uvalue: %f deg, p0: %f, q0: %f"%(optimized_uvalue,p0,q0))
    print("Radius: %f m "%radius)