#
# to profile it:
#

#
# pip install snakeviz
# pip install cprofilev
# pip install --upgrade pip
# python -B -m cProfile -o output.prof IR_WIG_shadow4_scan.py
# snakeviz output.prof

#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import numpy
import time
from srxraylib.plot.gol import plot,plot_scatter, set_qt

from shadow4.syned.magnetic_structure_1D_field import MagneticStructure1DField

from syned.storage_ring.electron_beam import ElectronBeam

from shadow4.sources.wiggler.source_wiggler import SourceWiggler
from shadow4.compatibility.beam3 import Beam3

from shadow4.beam.beam import Beam

import Shadow

from scipy.ndimage import gaussian_filter1d

def get_R(p_foc,q_foc,incidence):
    mm = (1.0 / p_foc + 1.0 / q_foc)
    return 2 / (numpy.cos(incidence * numpy.pi / 180)) / mm

def get_Rsagittal(p_foc,q_foc,incidence):
    mm = (1.0 / p_foc + 1.0 / q_foc)
    return 2 * (numpy.cos(incidence * numpy.pi / 180)) / mm

def run_source_wiggler_shadow4(electron_energy=2.0,
                               filename ="/home/manuel/Oasys/BM_only7.b",
                               use_emittances=True,
                               e_min = 1000.0,
                               e_max = 1000.1,
                               NRAYS = 5000, ):


    wigFile = "xshwig.sha"
    inData = ""

    nTrajPoints = 501
    ener_gev = electron_energy
    per = 0.5
    kValue = 4
    trajFile = ""
    shift_x_flag = 0
    shift_x_value = 0.0
    shift_betax_flag = 0
    shift_betax_value = 0.0


    sw = SourceWiggler()

    #
    # syned
    #


    syned_electron_beam = ElectronBeam(energy_in_GeV=electron_energy,
                                       current=0.5,
                                       moment_xx=   (7e-6          )**2    , #(39e-6)**2,
                                       moment_xpxp= (70e-12 / 7e-6 )**2    , #(2000e-12 / 51e-6)**2,
                                       moment_yy=   (10e-6         )**2    , #(31e-6)**2,
                                       moment_ypyp= (70e-12 / 10e-6)**2    , #(30e-12 / 31e-6)**2,
                                       )

    # conventional wiggler
    # syned_wiggler = Wiggler(K_vertical=kValue,K_horizontal=0.0,period_length=per,number_of_periods=nPer)

    # B from file


    syned_wiggler = MagneticStructure1DField.initialize_from_file(filename)
    # syned_wiggler.add_spatial_shift(-0.478)
    # syned_wiggler.flip_B()



    if e_min == e_max:
        ng_e = 1
    else:
        ng_e = 10

    sourcewiggler = SourceWiggler(name="test",
                    syned_electron_beam=syned_electron_beam,
                    syned_wiggler=syned_wiggler,
                    flag_emittance=use_emittances,
                    emin=e_min,
                    emax=e_max,
                    ng_e=ng_e,
                    ng_j=nTrajPoints)


    sourcewiggler.set_electron_initial_conditions_by_label(velocity_label="value_at_zero",
                                                           position_label="value_at_zero",)


    print(sourcewiggler.info())

    t00 = time.time()
    rays = sourcewiggler.calculate_rays(NRAYS=NRAYS)
    t11 = time.time() - t00
    print(">>>> time for %d rays: %f s, %f min, " % (NRAYS, t11, t11 / 60))

    beam = Beam3.initialize_from_array(rays)


    return beam




def run_source_wiggler_shadow3():
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



def run_beamline_moreno(beam, incidence=45.0, radius=1.0):

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


def run_beamline_toroid(beam,incidence=45.0,radius_major=10.0, radius_minor=1.0, p=10.0, q=10.0):

    iwrite = 0

    oe1 = Shadow.OE()
    oe2 = Shadow.OE()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #


    oe1.ALPHA = 90.0
    oe1.DUMMY = 100.0
    oe1.FMIRR = 3
    oe1.FWRITE = 1
    oe1.F_EXT = 1
    oe1.R_MAJ = radius_major
    oe1.R_MIN = radius_minor
    oe1.T_INCIDENCE = incidence
    oe1.T_REFLECTION = incidence
    oe1.T_SOURCE = p # 1.58 # 2.935
    oe1.T_IMAGE = q # 4.29 # 2.935

    oe2.ALPHA = 90.0
    oe2.DUMMY = 100.0
    oe2.FWRITE = 1
    oe2.F_REFRAC = 2
    oe2.T_IMAGE = 0.0
    oe2.T_INCIDENCE = 0.0
    oe2.T_REFLECTION = 180.0
    oe2.T_SOURCE = 0.0



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


def flux_vs_aperture(beam,nolost=1,rmax=None):
    x = beam.getshonecol(1,nolost=nolost)
    z = beam.getshonecol(2, nolost=nolost)
    r = numpy.sqrt(x**2 + z**2)
    intensity = beam.getshonecol(23, nolost=nolost)

    npoints = 100
    if rmax is None:
        rmax = 1.1*r.max()
    R = numpy.linspace(0,rmax,npoints)
    INT = numpy.zeros_like(R)
    for i in range(npoints):
        igood = numpy.argwhere( r <= R[i])
        if len(igood) == 0:
            INT[i] = 0
        else:
            INT[i] = intensity[igood].sum()
        print(">>>>>>", i, len(igood),INT[i],intensity.sum())

    plot(R,INT,xlog=1)






if __name__ == "__main__":

    from srxraylib.plot.gol import plot, set_qt
    from srxraylib.util.h5_simple_writer import H5SimpleWriter

    set_qt()


    do_calculate_source = False
    do_scan = False
    write_h5 = False
    do_plot = False
    all_magnets = False
    mirror_config = "toroid"  # moreno or toroid

    #
    #
    #


    if write_h5:
        h = H5SimpleWriter.initialize_file("IR_WIG_shadow4_scan.h5")



    if do_calculate_source:
        beam = run_source_wiggler_shadow4(electron_energy=2.0,
                                          filename ="/home/manuel/Oasys/BM_multi7.b",
                                          use_emittances=True,
                                          e_min = 0.01,
                                          e_max = 0.01,
                                          NRAYS = 100000, )

        beam.write("/home/manuel/Oasys/begin_ir_0p01_100000.dat")
    else:
        beam = Shadow.Beam()
        beam.load("/home/manuel/Oasys/begin_ir_0p01_100000.dat")


    if all_magnets == False:
        y = beam.rays[:, 1]

        ibad = numpy.where(y < -0.3)
        beam.rays[ibad, 6:9] = 0.0
        beam.rays[ibad, 15:18] = 0.0

        ibad = numpy.where(y > 0.3)
        beam.rays[ibad, 6:9] = 0.0
        beam.rays[ibad, 15:18] = 0.0

    beam_source = beam.duplicate()

    if do_plot:
        Shadow.ShadowTools.plotxy(beam,2,1,nbins=200)

    if mirror_config == "moreno":
        p_foc = 1.58
        q_foc = 4.290000
        Incidence = numpy.linspace(70, 80, 151)
    elif mirror_config == "toroid":
        p_foc =  2.935 # 5.87 / 2 # 1.58 #
        q_foc =  2.935 # 5.87 / 2 # 4.29 #
        # Incidence = numpy.linspace(73, 77, 151)
        # Incidence = numpy.linspace(70, 89, 151)
        Incidence = numpy.linspace(12, 88, 151)




    Radius   = numpy.zeros_like(Incidence)
    Fwhm     = numpy.zeros_like(Incidence)
    Std      = numpy.zeros_like(Incidence)
    Position = numpy.zeros_like(Incidence)

    if do_scan:
        for i,incidence in enumerate(Incidence):


            R = get_R(p_foc,q_foc,incidence)

            beam = None
            beam = beam_source.duplicate()

            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> R: ",R)

            if mirror_config == "moreno":
                beam, oe1 = run_beamline_moreno(beam, incidence=incidence, radius=R)
            elif mirror_config == "toroid":
                Rsag = get_Rsagittal(p_foc, q_foc, incidence)
                beam, oe1 = run_beamline_toroid(beam, incidence=incidence, radius_major=R, radius_minor=Rsag,
                      p=p_foc,q=q_foc)


            tkt = beam.histo1(1,ref=23,nbins=501,nolost=1,xrange=[-4000e-6,4000e-6])
            if tkt["fwhm"] is None:
                tkt["fwhm"] = 0.0
            imax = tkt["histogram"].argmax()
            Position[i] = 1e6 * tkt["bin_center"][imax]

            print(incidence,oe1.RMIRR,1e6*tkt["fwhm"])
            Radius[i] = R
            Fwhm[i] = 1e6*tkt["fwhm"]
            Std[i] = 1e6 * beam.get_standard_deviation(1,nolost=1,ref=23)

            if write_h5:
                h.create_entry("iteration incidence %f" % incidence, nx_default="histogram")
                h.add_dataset(1e6*tkt["bin_path"],tkt["histogram_path"], dataset_name="histogram",entry_name="iteration incidence %f" % incidence,
                              title_x="X / um",title_y="intensity / a.u.")



        if do_plot:
            plot(Incidence, Radius ,xtitle="Incidence/deg",ytitle="R/m",show=False)
            plot(Incidence, Fwhm,Incidence, Std,xtitle="Incidence/deg",ytitle="FWHM/um",legend=["fwhm","std"],show=False)
            plot(Incidence, Position, xtitle="Incidence/deg", ytitle="Position/um")

        if write_h5:
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
    else:
        optimizedIncidence = 45.0 # 85.75 # 74.946667 - 30e-2


    optimizedRadius    = get_R(p_foc=p_foc,q_foc=q_foc,incidence=optimizedIncidence)


    beam = None
    beam = beam_source.duplicate()
    if mirror_config == "moreno":
        beam, oe1 = run_beamline_moreno(beam, incidence=optimizedIncidence, radius=optimizedRadius)
    elif mirror_config == "toroid":
        Rsag = get_Rsagittal(p_foc, q_foc, optimizedIncidence)
        beam, oe1 = run_beamline_toroid(beam, incidence=optimizedIncidence, radius_major=optimizedRadius, radius_minor=Rsag,
                                        p=p_foc,q=q_foc)

    if do_plot:
        tkt = Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=201, nolost=1, ref=23,
                                        xrange=[-200e-6,200e-6],
                                        title="theta: %f, R: %f" % (optimizedIncidence, optimizedRadius))

    print("best incidence angle: %f deg, grazing %f deg" % (optimizedIncidence, 90 - optimizedIncidence))
    print("Radius: %f m " % optimizedRadius)

    if mirror_config == "toroid":

        print("Sagittal Radius: %f m " % Rsag)
        r1 = get_R(p_foc, q_foc, 45)
        r2 = get_Rsagittal(p_foc, q_foc, 45)
        print("Radii ar 45 deg: ",r1,r2)


    flux_vs_aperture(beam)