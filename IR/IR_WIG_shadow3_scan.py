#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow
import numpy


from twocylinders import shadow3file_twocylinders

def get_R(p_foc,q_foc,incidence):
    mm = (1.0 / p_foc + 1.0 / q_foc)
    return 2 / (numpy.cos(incidence * numpy.pi / 180)) / mm

def run_source_bm(NPOINT=30000,y_shift=0.0,magnetic_radius=10.0,horizontal_divergence=66e-3,
                  electron_energy=1.9):

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

    oe0.BENER = electron_energy
    if electron_energy == 1.9:
        oe0.EPSI_X = 1.989e-09
        oe0.EPSI_Z = 3.007e-11
        oe0.SIGMAX = 3.9e-05
        oe0.SIGMAZ = 3.1e-05
    else:
        oe0.EPSI_X = 70e-12
        oe0.EPSI_Z = 70e-12
        oe0.SIGMAX = 7e-06
        oe0.SIGMAZ = 10e-06
    oe0.FDISTR = 6
    oe0.FSOURCE_DEPTH = 4
    oe0.F_COLOR = 3
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.5 * horizontal_divergence
    oe0.HDIV2 = 0.5 * horizontal_divergence
    oe0.ISTAR1 = 5676561
    oe0.NCOL = 0
    oe0.NPOINT = NPOINT
    oe0.N_COLOR = 0
    oe0.PH1 = 0.4
    oe0.PH2 = 0.401
    oe0.POL_DEG = 0.0
    oe0.R_ALADDIN = magnetic_radius
    oe0.R_MAGNET = magnetic_radius
    oe0.SIGDIX = 0.0
    oe0.SIGDIZ = 0.0
    oe0.SIGMAY = 0.0
    oe0.VDIV1 = 0.05
    oe0.VDIV2 = 0.05
    oe0.WXSOU = 0.0
    oe0.WYSOU = 0.0
    oe0.WZSOU = 0.0


    beam.genSource(oe0)

    beam.rays[:, 1] += y_shift

    return beam

def run_adaptive_beamline(beam, incidence=45.0, FILE_RIP=b'/home/manuel/OASYS1.2/alsu-scripts/IR/presurface.dat', y_shift=0.0):

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
    oe1.FILE_RIP = FILE_RIP
    oe1.FWRITE = 1
    oe1.F_G_S = 2
    oe1.F_RIPPLE = 1
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
    oe2.SSOUR = 1.58 +1.18 - y_shift
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

    # Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space")
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")

    return beam,oe1



def run_beamline(beam, incidence=45.0, radius=1.0, y_shift=0.0):

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
    oe2.SSOUR = 1.58 +1.18 - y_shift
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

    # Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space")
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")

    return beam,oe1

def run_preprocessor(enerMin=1000.0,enerMax=1000.1,electron_energy=1.9):
    #
    # script to run the wiggler preprocessor (created by ShadowOui:Wiggler)
    #
    from srxraylib.sources import srfunc

    # (traj, pars) = srfunc.wiggler_trajectory(
    #     b_from=1,
    #     inData="BM_multi.b",
    #     nPer=1,
    #     nTrajPoints=501,
    #     ener_gev=1.9,
    #     per=0.01,
    #     kValue=1.0,
    #     trajFile="tmp.traj",
    #     shift_x_flag=4,
    #     shift_x_value=0.0,
    #     shift_betax_flag=4,
    #     shift_betax_value=0.094)

    (traj, pars) = srfunc.wiggler_trajectory(
        b_from=1,
        inData="BM_multi.b",
        nPer=1,
        nTrajPoints=501,
        ener_gev=electron_energy,
        per=0.01,
        kValue=1.0,
        trajFile="tmp.traj",
        shift_x_flag=5,
        shift_x_value=0.042,
        shift_betax_flag=4,
        shift_betax_value=0.0324)

    #
    # calculate cdf and write file for Shadow/Source
    #

    srfunc.wiggler_cdf(traj,
                       enerMin=enerMin,
                       enerMax=enerMax,
                       enerPoints=1001,
                       outFile=b'xshwig.sha',
                       elliptical=False)

    calculate_spectrum = False

    if calculate_spectrum:
        e, f, w = srfunc.wiggler_spectrum(traj,
                                          enerMin=1000.0,
                                          enerMax=100000.1,
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



def run_source_wiggler(select_bm=1,electron_energy=1.9):
    #
    # initialize shadow3 source (oe0) and beam
    #

    oe0 = Shadow.Source()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.BENER = electron_energy
    if electron_energy == 1.9:
        oe0.EPSI_X = 1.989e-09
        oe0.EPSI_Z = 3.007e-11
        oe0.SIGMAX = 3.9e-05
        oe0.SIGMAZ = 3.1e-05
    else:
        oe0.EPSI_X = 70e-12
        oe0.EPSI_Z = 70e-12
        oe0.SIGMAX = 7e-06
        oe0.SIGMAZ = 10e-06
    oe0.CONV_FACT = 1.0
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
    oe0.SIGMAY = 0.0
    oe0.VDIV1 = 1.0
    oe0.VDIV2 = 1.0
    oe0.WXSOU = 0.0
    oe0.WYSOU = 0.0
    oe0.WZSOU = 0.0


    # Run SHADOW to create the source

    beam = Shadow.Beam()

    beam.genSource(oe0)

    #
    # select a part
    #
    if select_bm == 1:
        y = beam.rays[:, 1]
        ibad = numpy.where(y > -0.2)
        beam.rays[ibad, 6:9] = 0.0
        beam.rays[ibad, 15:18] = 0.0
    elif select_bm == 2:
        pass
    elif select_bm == 3:
        y = beam.rays[:, 1]
        ibad = numpy.where(y < 0.2)
        beam.rays[ibad, 6:9] = 0.0
        beam.rays[ibad, 15:18] = 0.0

    return beam

if __name__ == "__main__":
    import h5py
    from srxraylib.plot.gol import plot, set_qt
    from srxraylib.util.h5_simple_writer import H5SimpleWriter

    set_qt()

    wiggler_or_bm = 4 # 0=wiggler, 1=Mag7, 2=Antibend, 3=Mag8, 4=ALS
    select_bm = 1 # this is only for wiggler selection: 1=Mag7 2=antibend(to check), 3=Mag8
    use_adaptive = True


    Incidence = numpy.linspace(55, 80, 151)
    # Incidence = numpy.linspace(80, 89, 50)

    Radius   = numpy.zeros_like(Incidence)
    Fwhm     = numpy.zeros_like(Incidence)
    Std      = numpy.zeros_like(Incidence)
    Position = numpy.zeros_like(Incidence)


    h = H5SimpleWriter.initialize_file("IR_WIG_shadow3_scan.h5")

    # Radius M1: 7.615618611829955
    # Radius AB: 41.695511899769
    # Radius M2: 7.857811429874018
    # Half - Divergence M1: 0.032827274151
    # Half - Divergence AB: 0.0038973019540000007
    # Half - Divergence M2: 0.031841706445325

    if wiggler_or_bm == 0:
        run_preprocessor(enerMin=1000.0, enerMax=1001.0)
        beam = run_source_wiggler(select_bm=select_bm,electron_energy=2.0)

        y = beam.getshonecol(2)
        w = beam.getshonecol(23)
        y_shift = numpy.average(y, weights=w)
    elif wiggler_or_bm == 1:
        y_shift =  0.0 # -0.4778 # 0.0
        magnetic_radius=-7.615618611829955
        horizontal_divergence = 2*0.032827274151
        beam = run_source_bm(y_shift=y_shift,magnetic_radius=magnetic_radius,
                             horizontal_divergence=horizontal_divergence,
                             electron_energy=2.0)
    elif wiggler_or_bm == 2:
        y_shift =  0.0
        magnetic_radius= 41.695511899769 ## attention to sign reversed!!
        horizontal_divergence = 2*0.0038973019540000007
        beam = run_source_bm(y_shift=y_shift,magnetic_radius=magnetic_radius,
                             horizontal_divergence=horizontal_divergence,
                             electron_energy=2.0)
    elif wiggler_or_bm == 3:
        y_shift =  0.0 # 0.4778 # 0.0
        magnetic_radius=-7.857811429874018
        horizontal_divergence = 2*0.031841706445325
        beam = run_source_bm(y_shift=y_shift,magnetic_radius=magnetic_radius,
                             horizontal_divergence=horizontal_divergence,
                             electron_energy=2.0)
    elif wiggler_or_bm == 4:
        y_shift = 0.0  # 0.0
        magnetic_radius = -5.0
        horizontal_divergence = 69e-3
        beam = run_source_bm(y_shift=y_shift, magnetic_radius=magnetic_radius,
                             horizontal_divergence=horizontal_divergence,
                             electron_energy=1.9)

    print(">>>>>y_shift: ",y_shift)

    beam_source = beam.duplicate()


    for i,incidence in enumerate(Incidence):

        p_foc = 1.58 - y_shift  # 2.058216
        q_foc = 4.290000
        # mm = (1.0 / p_foc + 1.0 / q_foc)
        # R = 2 / (numpy.cos(incidence * numpy.pi / 180)) / mm
        R = get_R(p_foc,q_foc,incidence)

        # if bm_or_wiggler == 0:
        beam = None
        beam = beam_source.duplicate()

        R1 = R # get_R(p_foc + 0.5 , q_foc, incidence)
        R2 = R # get_R(p_foc - 0.5 , q_foc, incidence)

        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> R: ",R)
        if use_adaptive:
            shadow3file_twocylinders(numpy.linspace(-0.02, 0.02, 11),
                                     numpy.linspace(-0.25, 0.25, 200),
                                     delta=0.00, radius1=R1, radius2=R2, filename="presurface.dat")
            beam, oe1 = run_adaptive_beamline(beam, incidence=incidence, FILE_RIP=b"presurface.dat", y_shift=y_shift)
        else:
            beam, oe1 = run_beamline(beam, incidence=incidence, radius=R, y_shift=y_shift)



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

    #
    # plot best result
    #
    if wiggler_or_bm == 0:
        beam = None
        beam = beam_source.duplicate()
    else:
        beam = run_source_bm(y_shift=y_shift, magnetic_radius=magnetic_radius,
                             horizontal_divergence=horizontal_divergence,NPOINT=300000)
        # beam = run_source_bm(y_shift=y_shift, magnetic_radius=magnetic_radius, NPOINT=300000)

    if use_adaptive:
        shadow3file_twocylinders(numpy.linspace(-0.02, 0.02, 11),
                                 numpy.linspace(-0.25, 0.25, 200),
                                 delta=0.00, radius1=optimizedRadius, radius2=optimizedRadius, filename="presurface.dat")
        beam, oe1 = run_adaptive_beamline(beam, incidence=optimizedIncidence, FILE_RIP=b"presurface.dat", y_shift=y_shift)
    else:
        beam, oe1 = run_beamline(beam, incidence=optimizedIncidence, radius=optimizedRadius, y_shift=y_shift)


    print("best incidence angle: %f deg, grazing: %f deg"%(optimizedIncidence,90-optimizedIncidence))
    print("Radius: %f m"%optimizedRadius)
    tkt = Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=201, nolost=1, ref=23,
                              # xrange=[-5e-3*5,5e-3*5],
                              title="theta: %f, R: %f"%(optimizedIncidence,optimizedRadius))

    for k in tkt.keys():
        print(">>> k: ",k, tkt[k])

    print("y_shift: ", y_shift)
    print("best incidence angle: %f deg, grazing: %f deg"%(optimizedIncidence,90-optimizedIncidence))
    print("Radius: %f m"%optimizedRadius)
    try:
        print("Focal size FWHM: H: %f um, V: %f um: "%(1e6*tkt["fwhm_h"],1e6*tkt["fwhm_v"]))
    except:
        pass


    # #
    # # optimize shift
    # #
    # Shift = numpy.linspace(0.9,1.1,100)
    # Radius = optimizedRadius * Shift #numpy.zeros_like(Shift)
    # Fwhm = numpy.zeros_like(Shift)
    # Std = numpy.zeros_like(Shift)
    # SPosition = numpy.zeros_like(Shift)
    # Position = numpy.zeros_like(Shift)
    #
    # incidence = optimizedIncidence
    #
    # for i in range(Shift.size):
    #
    #     beam = None
    #     beam = beam_source.duplicate()
    #
    #     beam0, oe1 = run_shadow(source=False, trace=True, p_shift=p_shift, beam=beam,
    #                             incidence=incidence,radius=Radius[i])
    #
    #
    #     # Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space")
    #
    #     tkt = beam.histo1(1,ref=23,nbins=201,nolost=1) # xrange=[-4000e-6,4000e-6],
    #     if tkt["fwhm"] is None:
    #         tkt["fwhm"] = 0.0
    #     imax = tkt["histogram"].argmax()
    #
    #     print(incidence,oe1.RMIRR,1e6*tkt["fwhm"])
    #     # Radius[i] = oe1.RMIRR
    #     Fwhm[i] = 1e6*tkt["fwhm"]
    #     Std[i] = 1e6 * beam.get_standard_deviation(1,nolost=1,ref=23)
    #
    #     tkt = beam.histo1(1,ref=23,nolost=1,xrange=[-5e-3,5e-3])
    #     imax = tkt["histogram"].argmax()
    #     Position[i] = tkt["bin_center"][imax]
    #     #
    #     # h.create_entry("iteration incidence %f" % incidence, nx_default="histogram")
    #     # h.add_dataset(1e6*tkt["bin_path"],tkt["histogram_path"], dataset_name="histogram",entry_name="iteration incidence %f" % incidence,
    #     #               title_x="X / um",title_y="intensity / a.u.")
    #     #
    # plot(Radius, Fwhm,
    #      Radius, Std,xtitle="radius/m",ytitle="FWHM/um",legend=["fwhm","std"],show=False)
    # plot(Radius, Position,       xtitle="radius/m", ytitle="Position/um")
    #
    # print("best Radius before optimization: ", optimizedRadius)
    # imin = Std.argmin()
    # optimizedRadius = Radius[imin]
    # print("best Radius: ", optimizedRadius)