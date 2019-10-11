#
# Python script to run shadow4. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow
import numpy


def run_shadow(incidence=67.7,source=True,trace=True,beam=None):
    #
    # Python script to run shadow4. Created automatically with ShadowTools.make_python_script_from_list().
    #
    import Shadow
    import numpy
    from shadow4.compatibility.beam3 import Beam3

    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 0

    #
    # initialize shadow4 source (oe0) and beam
    #
    # beam = Beam3()
    oe0 = Shadow.Source()
    oe1 = Shadow.OE()
    oe2 = Shadow.OE()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.BENER = 1.9
    oe0.EPSI_X = 1.989e-09
    oe0.EPSI_Z = 3.007e-11
    oe0.FDISTR = 6
    oe0.FSOURCE_DEPTH = 4
    oe0.F_COLOR = 3
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.035
    oe0.HDIV2 = 0.035
    oe0.ISTAR1 = 5676561
    oe0.NCOL = 0
    oe0.NPOINT = 50000
    oe0.N_COLOR = 0
    oe0.PH1 = 0.4
    oe0.PH2 = 0.401
    oe0.POL_DEG = 0.0
    oe0.R_ALADDIN = -5.0
    oe0.R_MAGNET = -5.0
    oe0.SIGDIX = 0.0
    oe0.SIGDIZ = 0.0
    oe0.SIGMAX = 3.9e-05
    oe0.SIGMAY = 0.0
    oe0.SIGMAZ = 3.1e-05
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
    oe1.F_DEFAULT = 0
    oe1.SIMAG = 4.29
    oe1.SSOUR = 1.58
    oe1.THETA = incidence
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

    # Run SHADOW to create the source
    if source:
        if iwrite:
            oe0.write("start.00")

        # beam.genSource(oe0)
        from shadow4.compatibility.translate_start_files import start00_to_bending_magnet
        bm4 = start00_to_bending_magnet(oe0)
        print(">>>>calling calculate_rays...")
        t00 = time.time()
        rays = bm4.calculate_rays(F_COHER=oe0.F_COHER,NRAYS=oe0.NPOINT,SEED=oe0.ISTAR1,
                                     EPSI_DX=oe0.EPSI_DX,EPSI_DZ=oe0.EPSI_DZ,verbose=False)
        t11 = time.time()-t00
        print(">>>> time for %d rays: %f s, %f min, "%(oe0.NPOINT,t11,t11/60))
        beam = Beam3.initialize_from_array(rays)

        if iwrite:
            oe0.write("end.00")
            beam.write("begin.dat")

    if not trace:
        return beam
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

if __name__ == "__main__":
    import time
    from srxraylib.plot.gol import plot, set_qt
    from srxraylib.util.h5_simple_writer import H5SimpleWriter

    t0 = time.time()
    set_qt()

    Grazing = numpy.linspace(15, 30.0, 151)

    Radius = numpy.zeros_like(Grazing)
    Fwhm = numpy.zeros_like(Grazing)
    Std = numpy.zeros_like(Grazing)

    import h5py
    h = H5SimpleWriter.initialize_file("IR_BM_shadow4_scan.h5")
    Incidence = 90.0 - Grazing

    beam = run_shadow(source=True,trace=False,beam=None)

    beam_source = beam.duplicate()


    for i,incidence in enumerate(Incidence):
        grazing = 90.0 - incidence



        # beam,oe1 = run_shadow(incidence=incidence)
        beam = None
        beam = beam_source.duplicate()
        beam,oe1 = run_shadow(incidence=incidence, source=False, trace=True, beam=beam)


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

    # beam = run_shadow(source=True,trace=False,beam=None)
    #
    # beam_source = beam.duplicate()
    #
    # beam = None
    # beam = beam_source.duplicate()
    # beam,oe1 = run_shadow(incidence=30.0,source=False, trace=True, beam=beam)
    # Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space %f"%oe1.T_INCIDENCE)
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")
    #
    #
    # beam = None
    # beam = beam_source.duplicate()
    # beam,oe1 = run_shadow(incidence=67.7,source=False, trace=True, beam=beam)
    # Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space %f"%oe1.T_INCIDENCE)
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")


    print(">>Time spent: ",time.time()-t0)