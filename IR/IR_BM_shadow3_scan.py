#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow
import numpy

def get_R(p_foc,q_foc,incidence):
    mm = (1.0 / p_foc + 1.0 / q_foc)
    return 2 / (numpy.cos(incidence * numpy.pi / 180)) / mm


def run_shadow(incidence=67.7,radius=1.0):
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

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.BENER = 2.0
    oe0.EPSI_X = 7e-11
    oe0.EPSI_Z = 7e-11
    oe0.FDISTR = 6
    oe0.FSOURCE_DEPTH = 4
    oe0.F_COLOR = 3
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.033
    oe0.HDIV2 = 0.033
    oe0.ISTAR1 = 5676561
    oe0.NCOL = 0
    oe0.NPOINT = 300000
    oe0.N_COLOR = 0
    oe0.PH1 = 0.1
    oe0.PH2 = 0.101
    oe0.POL_DEG = 0.0
    oe0.R_ALADDIN = -7.615618611829955
    oe0.R_MAGNET = -7.615618611829955
    oe0.SIGDIX = 0.0
    oe0.SIGDIZ = 0.0
    oe0.SIGMAX = 7e-06
    oe0.SIGMAY = 0.0
    oe0.SIGMAZ = 1e-05
    oe0.VDIV1 = 0.05
    oe0.VDIV2 = 0.05
    oe0.WXSOU = 0.0
    oe0.WYSOU = 0.0
    oe0.WZSOU = 0.0

    oe1.ALPHA = 90.0
    oe1.DUMMY = 100.0
    oe1.FCYL = 1
    oe1.FMIRR = 1
    oe1.FWRITE = 1
    oe1.F_EXT = 1
    oe1.RMIRR = radius
    oe1.T_IMAGE = 0.0
    oe1.T_INCIDENCE = incidence
    oe1.T_REFLECTION = incidence
    oe1.T_SOURCE = 1.58

    oe2.ALPHA = 90.0
    oe2.DUMMY = 100.0
    oe2.FCYL = 1
    oe2.FHIT_C = 1
    oe2.FMIRR = 1
    oe2.F_DEFAULT = 0
    oe2.RLEN1 = 6.0
    oe2.RLEN2 = 6.0
    oe2.RWIDX1 = 0.65
    oe2.RWIDX2 = 0.65
    oe2.SIMAG = 2.697
    oe2.SSOUR = 3.31
    oe2.THETA = 45.0
    oe2.T_IMAGE = 2.697
    oe2.T_INCIDENCE = 45.0
    oe2.T_REFLECTION = 45.0
    oe2.T_SOURCE = 1.73

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

    # Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space")
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")

    return beam, oe1

if __name__ == "__main__":

    from srxraylib.plot.gol import plot, set_qt
    from srxraylib.util.h5_simple_writer import H5SimpleWriter

    set_qt()

    Grazing = numpy.linspace(10, 20.0, 151)

    Radius = numpy.zeros_like(Grazing)
    Fwhm = numpy.zeros_like(Grazing)
    Std = numpy.zeros_like(Grazing)

    import h5py

    h = H5SimpleWriter.initialize_file("IR_BM_shadow3_scan.h5")

    Incidence = 90.0 - Grazing
    
    for i,incidence in enumerate(Incidence):
        grazing = 90.0 - incidence
        radius = get_R(1.58,6.007-1.58,incidence)
        beam,oe1 = run_shadow(incidence=incidence,radius=radius)
        # Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space")

        tkt = beam.histo1(1,nbins=201,nolost=1) # xrange=[-4000e-6,4000e-6],
        print(incidence,oe1.RMIRR,1e6*tkt["fwhm"])
        Radius[i] = oe1.RMIRR
        Fwhm[i] = 1e6*tkt["fwhm"]
        Std[i] = 1e6 * beam.get_standard_deviation(1,nolost=1)

        h.create_entry("iteration %f" % grazing, nx_default="histogram")
        h.add_dataset(1e6*tkt["bin_path"],tkt["histogram_path"], dataset_name="histogram",entry_name="iteration %f" % grazing,
                      title_x="X / um",title_y="intensity / a.u.")



    plot(Radius,Fwhm,xtitle="R/m",ytitle="FWHM/um",show=False)
    plot(Grazing, Fwhm,Grazing, Std,xtitle="Graz/deg",ytitle="FWHM/um",legend=["fwhm","std"])

    h.create_entry("scan results", nx_default="Fwhm")
    h.add_dataset(Grazing, Fwhm, dataset_name="Fwhm", entry_name="scan results",
                  title_x="grazing angle / deg", title_y="fwhm / um")
    h.add_dataset(Grazing, Std, dataset_name="Std", entry_name="scan results",
                  title_x="grazing angle / deg", title_y="std / um")
    h.add_dataset(Radius, Fwhm, dataset_name="Radius", entry_name="scan results",
                  title_x="radius / m", title_y="fwhm / um")

    for key in tkt.keys():
        print(key)
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")
    imin = Std.argmin()
    optimizedIncidence = Incidence[imin]     #
    optimizedRadius    = Radius[imin]        #



    beam,oe1 = run_shadow(incidence=optimizedIncidence,radius=optimizedRadius)


    tkt = Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=201, nolost=1, ref=23,
                              # xrange=[-5e-3*5,5e-3*5],
                              title="theta: %f, R: %f"%(optimizedIncidence,optimizedRadius))

    print("best incidence angle: %f deg, grazing %f deg"%(optimizedIncidence,90-optimizedIncidence))
    print("Radius: %f m "%optimizedRadius)