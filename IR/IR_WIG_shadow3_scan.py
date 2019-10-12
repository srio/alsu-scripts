#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow
import numpy


def run_preprocessor():
    #
    # script to run the wiggler preprocessor (created by ShadowOui:Wiggler)
    #
    from srxraylib.sources import srfunc

    (traj, pars) = srfunc.wiggler_trajectory(
        b_from=1,
        inData="BM_multi.b",
        nPer=1,
        nTrajPoints=501,
        ener_gev=1.9,
        per=0.01,
        kValue=1.0,
        trajFile="tmp.traj",
        shift_x_flag=4,
        shift_x_value=0.0,
        shift_betax_flag=4,
        shift_betax_value=0.094)

    #
    # calculate cdf and write file for Shadow/Source
    #

    srfunc.wiggler_cdf(traj,
                       enerMin=1000.0,
                       enerMax=1000.1,
                       enerPoints=1001,
                       outFile=b'xshwig.sha',
                       elliptical=False)

    calculate_spectrum = False

    if calculate_spectrum:
        e, f, w = srfunc.wiggler_spectrum(traj,
                                          enerMin=1000.0,
                                          enerMax=1000.1,
                                          nPoints=500,
                                          electronCurrent=400.0 * 1e-3,
                                          outFile="spectrum.dat",
                                          elliptical=False)
        from srxraylib.plot.gol import plot
        plot(e, f, xlog=False, ylog=False, show=False,
             xtitle="Photon energy [eV]", ytitle="Flux [Photons/s/0.1%bw]", title="Flux")
        plot(e, w, xlog=False, ylog=False, show=True,
             xtitle="Photon energy [eV]", ytitle="Spectral Power [E/eV]", title="Spectral Power")
        #
    # end script
    #


def run_shadow(source=True, obscure=3, trace=True, beam=None, p_shift=0.0, incidence=39.9):
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

    oe0 = Shadow.Source()
    oe1 = Shadow.OE()
    oe2 = Shadow.OE()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.BENER = 1.9
    oe0.CONV_FACT = 1.0
    oe0.EPSI_X = 1.989e-09
    oe0.EPSI_Z = 3.007e-11
    oe0.FDISTR = 0
    oe0.FILE_TRAJ = b'xshwig.sha'
    oe0.FSOUR = 0
    oe0.FSOURCE_DEPTH = 0
    oe0.F_COLOR = 0
    oe0.F_PHOT = 0
    oe0.F_WIGGLER = 1
    oe0.HDIV1 = 1.0
    oe0.HDIV2 = 1.0
    oe0.IDO_VX = 0
    oe0.IDO_VZ = 0
    oe0.IDO_X_S = 0
    oe0.IDO_Y_S = 0
    oe0.IDO_Z_S = 0
    oe0.ISTAR1 = 5676561
    oe0.NCOL = 0
    oe0.NPOINT = 50000
    oe0.N_COLOR = 0
    oe0.PH1 = 1000.0
    oe0.PH2 = 1000.1
    oe0.POL_DEG = 0.0
    oe0.SIGMAX = 3.9e-05
    oe0.SIGMAY = 0.0
    oe0.SIGMAZ = 3.1e-05
    oe0.VDIV1 = 1.0
    oe0.VDIV2 = 1.0
    oe0.WXSOU = 0.0
    oe0.WYSOU = 0.0
    oe0.WZSOU = 0.0

    oe1.ALPHA = 90.0
    oe1.DUMMY = 100.0
    oe1.FCYL = 1
    oe1.FHIT_C = 1
    oe1.FMIRR = 1
    oe1.FWRITE = 1
    oe1.F_DEFAULT = 0
    oe1.RLEN1 = 0.3
    oe1.RLEN2 = 0.4
    oe1.RWIDX1 = 0.65
    oe1.RWIDX2 = 0.65
    oe1.SIMAG = 4.29
    oe1.THETA = 39.0
    oe1.T_IMAGE = 0.0
    oe1.T_INCIDENCE = incidence
    oe1.T_REFLECTION = incidence
    oe1.T_SOURCE = 1.58
    oe1.SSOUR = oe1.T_SOURCE - p_shift

    oe2.ALPHA = 90.0
    oe2.DUMMY = 100.0
    oe2.FCYL = 1
    oe2.FMIRR = 1
    oe2.FWRITE = 1
    oe2.F_DEFAULT = 0
    oe2.SIMAG = 3.11
    oe2.SSOUR = 2.76
    oe2.THETA = 45.0
    oe2.T_IMAGE = 3.11
    oe2.T_INCIDENCE = 45.0
    oe2.T_REFLECTION = 45.0
    oe2.T_SOURCE = 1.18

    # Run SHADOW to create the source

    if source:
        beam = Shadow.Beam()
        if iwrite:
            oe0.write("start.00")

        beam.genSource(oe0)

        #
        # obscure a part
        #
        if obscure == 3:
            y = beam.rays[:, 1]
            ibad = numpy.where(y > -0.2)
            beam.rays[ibad, 6:9] = 0.0
            beam.rays[ibad, 15:18] = 0.0

        if iwrite:
            oe0.write("end.00")
            beam.write("begin.dat")

    if not trace:
        return beam, oe1




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

    Incidence = numpy.linspace(30, 50.0, 51)
    Grazing = 90.0 - Incidence

    Radius = numpy.zeros_like(Grazing)
    Fwhm = numpy.zeros_like(Grazing)
    Std = numpy.zeros_like(Grazing)
    Position = numpy.zeros_like(Grazing)

    import h5py
    print(">>>>>>>>>>>>>>>>>>>>",h5py.version.version)
    h = H5SimpleWriter.initialize_file("IR_WIG_shadow3_scan.h5")



    run_preprocessor()
    beam, oe1 = run_shadow(source=True,obscure=3,trace=False)
    beam_source = beam.duplicate()

    y = beam_source.getshonecol(2)
    w = beam_source.getshonecol(23)
    # tkt = beam.histo1(2, ref=23, nolost=1, nbins=301)
    # imax = tkt["histogram"].argmax()
    # center = tkt["bin_center"][imax]

    p_shift = numpy.average(y,weights=w)
    print(">>>>>o_shift: ",p_shift)

    # Shadow.ShadowTools.plotxy(beam, 2, 1, nbins=101, nolost=1, title="Top view")
    # Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space")
    
    for i,incidence in enumerate(Incidence):
        grazing = 90.0 - incidence
        beam = None
        beam = beam_source.duplicate()

        beam,oe1 = run_shadow(source=False, trace=True, p_shift=p_shift, beam=beam, incidence=incidence)
        # Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space")

        tkt = beam.histo1(1,ref=23,nbins=201,nolost=1) # xrange=[-4000e-6,4000e-6],
        if tkt["fwhm"] is None:
            tkt["fwhm"] = 0.0
        imax = tkt["histogram"].argmax()
        Position[i] = 1e6 * tkt["bin_center"][imax]

        print(incidence,oe1.RMIRR,1e6*tkt["fwhm"])
        Radius[i] = oe1.RMIRR
        Fwhm[i] = 1e6*tkt["fwhm"]
        Std[i] = 1e6 * beam.get_standard_deviation(1,nolost=1,ref=23)
        # tkt = beam.histo1(1,ref=23,nolost=1,xrange=[-5e-3,5e-3])
        # imax = tkt["histogram"].argmax()
        # Position[i] = tkt["bin_center"][imax]

        h.create_entry("iteration incidence %f" % incidence, nx_default="histogram")
        h.add_dataset(1e6*tkt["bin_path"],tkt["histogram_path"], dataset_name="histogram",entry_name="iteration incidence %f" % incidence,
                      title_x="X / um",title_y="intensity / a.u.")



    # plot(Radius,Fwhm,xtitle="R/m",ytitle="FWHM/um",show=False)
    # plot(Incidence, Fwhm,Incidence, Std,xtitle="Incidence/deg",ytitle="FWHM/um",legend=["fwhm","std"],show=False)
    # plot(Incidence, Position, xtitle="Incidence/deg", ytitle="Position/um")

    h.create_entry("scan results", nx_default="Fwhm")
    h.add_dataset(Incidence, Fwhm, dataset_name="Fwhm", entry_name="scan results",
                  title_x="incidence angle / deg", title_y="fwhm / um")
    h.add_dataset(Incidence, Std, dataset_name="Std", entry_name="scan results",
                  title_x="incidence angle / deg", title_y="std / um")
    h.add_dataset(Radius, Fwhm, dataset_name="Radius", entry_name="scan results",
                  title_x="radius / m", title_y="fwhm / um")

    h.add_dataset(Incidence, Position, dataset_name="PositionCenter", entry_name="scan results",
                  title_x="incidence angle / deg", title_y="center / um")

    for key in tkt.keys():
        print(key)

    print("p_shift: ",p_shift)
    print("p: ", oe1.SSOUR)
    imin = Std.argmin()
    print("best incidence angle: ",Incidence[imin])
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")

    #
    # optimize shift
    #
    Shift = numpy.linspace(0.8,1.2,100)
    Radius = numpy.zeros_like(Shift)
    Fwhm = numpy.zeros_like(Shift)
    Std = numpy.zeros_like(Shift)
    SPosition = numpy.zeros_like(Shift)
    Position = numpy.zeros_like(Shift)

    incidence = Incidence[imin]
    for i in range(Shift.size):
        beam = None
        beam = beam_source.duplicate()

        beam,oe1 = run_shadow(source=False, trace=True, p_shift=p_shift*Shift[i], beam=beam, incidence=incidence)
        # Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space")

        tkt = beam.histo1(1,ref=23,nbins=201,nolost=1) # xrange=[-4000e-6,4000e-6],
        if tkt["fwhm"] is None:
            tkt["fwhm"] = 0.0
        imax = tkt["histogram"].argmax()
        Position[i] = 1e6 * tkt["bin_center"][imax]

        print(incidence,oe1.RMIRR,1e6*tkt["fwhm"])
        Radius[i] = oe1.RMIRR
        Fwhm[i] = 1e6*tkt["fwhm"]
        Std[i] = 1e6 * beam.get_standard_deviation(1,nolost=1,ref=23)
        SPosition[i] = oe1.SSOUR
        # # tkt = beam.histo1(1,ref=23,nolost=1,xrange=[-5e-3,5e-3])
        # # imax = tkt["histogram"].argmax()
        # # Position[i] = tkt["bin_center"][imax]
        #
        # h.create_entry("iteration incidence %f" % incidence, nx_default="histogram")
        # h.add_dataset(1e6*tkt["bin_path"],tkt["histogram_path"], dataset_name="histogram",entry_name="iteration incidence %f" % incidence,
        #               title_x="X / um",title_y="intensity / a.u.")
        #
    plot(SPosition, Fwhm,xtitle="p/m",ytitle="FWHM/um",show=False)
    plot(SPosition, Fwhm,SPosition, Std,xtitle="p/m",ytitle="FWHM/um",legend=["fwhm","std"],show=False)
    plot(SPosition, Position, xtitle="p/m", ytitle="Position/um")
