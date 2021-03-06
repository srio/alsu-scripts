#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow
import numpy


def get_R(p_foc,q_foc,incidence):
    mm = (1.0 / p_foc + 1.0 / q_foc)
    return 2 / (numpy.cos(incidence * numpy.pi / 180)) / mm


def run_source_wiggler():
    from srxraylib.sources import srfunc
    (traj, pars) = srfunc.wiggler_trajectory(
        b_from=1,
        inData="/home/manuel/Oasys/BM_only7.b",
        nPer=1,
        nTrajPoints=501,
        ener_gev=2.0,
        per=0.01,
        kValue=1.0,
        trajFile="tmp.traj",
        shift_x_flag=4,
        shift_x_value=0.042,
        shift_betax_flag=4,
        shift_betax_value=0.035)

    #
    # calculate cdf and write file for Shadow/Source
    #

    srfunc.wiggler_cdf(traj,
                       enerMin=1000.0,
                       enerMax=1000.1,
                       enerPoints=1001,
                       outFile=b'/home/manuel/Oasys/xshwig.sha',
                       elliptical=False)

    calculate_spectrum = False

    if calculate_spectrum:
        e, f, w = srfunc.wiggler_spectrum(traj,
                                          enerMin=1000.0,
                                          enerMax=1000.1,
                                          nPoints=500,
                                          electronCurrent=500.0 * 1e-3,
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


    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 0

    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = Shadow.Beam()
    oe0 = Shadow.Source()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.BENER = 2.0
    oe0.CONV_FACT = 1.0
    oe0.EPSI_X = 70e-12
    oe0.EPSI_Z = 70e-12
    oe0.FDISTR = 0
    oe0.FILE_TRAJ = b'/home/manuel/Oasys/xshwig.sha'
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
    oe0.SIGMAX = 7e-06
    oe0.SIGMAY = 0.0
    oe0.SIGMAZ = 1e-05
    oe0.VDIV1 = 1.0
    oe0.VDIV2 = 1.0
    oe0.WXSOU = 0.0
    oe0.WYSOU = 0.0
    oe0.WZSOU = 0.0



    #Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")

    return beam



def run_beamline(beam, incidence=45.0, radius=1.0):

    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 0

    #
    # initialize shadow3 source (oe0) and beam
    #

    oe1 = Shadow.OE()
    oe2 = Shadow.OE()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

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

    return beam,oe1



if __name__ == "__main__":
    import h5py
    from srxraylib.plot.gol import plot, set_qt
    from srxraylib.util.h5_simple_writer import H5SimpleWriter

    set_qt()


    use_adaptive = True


    Incidence = numpy.linspace(70, 80, 151)


    Radius   = numpy.zeros_like(Incidence)
    Fwhm     = numpy.zeros_like(Incidence)
    Std      = numpy.zeros_like(Incidence)
    Position = numpy.zeros_like(Incidence)


    h = H5SimpleWriter.initialize_file("IR_WIG_shadow3_scan.h5")


    beam = run_source_wiggler()

    beam_source = beam.duplicate()

    Shadow.ShadowTools.plotxy(beam,2,1,nbins=200)


    for i,incidence in enumerate(Incidence):

        p_foc = 1.58
        q_foc = 4.290000
        R = get_R(p_foc,q_foc,incidence)

        beam = None
        beam = beam_source.duplicate()

        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> R: ",R)

        beam, oe1 = run_beamline(beam, incidence=incidence, radius=R)



        tkt = beam.histo1(1,ref=23,nbins=201,nolost=1) # xrange=[-4000e-6,4000e-6],
        if tkt["fwhm"] is None:
            tkt["fwhm"] = 0.0
        imax = tkt["histogram"].argmax()
        Position[i] = 1e6 * tkt["bin_center"][imax]

        print(incidence,oe1.RMIRR,1e6*tkt["fwhm"])
        Radius[i] = R
        Fwhm[i] = 1e6*tkt["fwhm"]
        Std[i] = 1e6 * beam.get_standard_deviation(1,nolost=1,ref=23)

        h.create_entry("iteration incidence %f" % incidence, nx_default="histogram")
        h.add_dataset(1e6*tkt["bin_path"],tkt["histogram_path"], dataset_name="histogram",entry_name="iteration incidence %f" % incidence,
                      title_x="X / um",title_y="intensity / a.u.")



    plot(Incidence, Radius ,xtitle="Incidence/deg",ytitle="R/m",show=False)
    plot(Incidence, Fwhm,Incidence, Std,xtitle="Incidence/deg",ytitle="FWHM/um",legend=["fwhm","std"],show=False)
    plot(Incidence, Position, xtitle="Incidence/deg", ytitle="Position/um")

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


    imin = Std.argmin()
    optimizedIncidence = Incidence[imin]     #
    optimizedRadius    = Radius[imin]        #

    beam = None
    beam = beam_source.duplicate()
    beam, oe1 = run_beamline(beam, incidence=optimizedIncidence, radius=optimizedRadius)

    tkt = Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=201, nolost=1, ref=23,
                                    # xrange=[-5e-3*5,5e-3*5],
                                    title="theta: %f, R: %f" % (optimizedIncidence, optimizedRadius))

    print("best incidence angle: %f deg, grazing %f deg" % (optimizedIncidence, 90 - optimizedIncidence))
    print("Radius: %f m " % optimizedRadius)
