#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow

def run_shadow3_source():
    #
    # Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
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
    oe0.NPOINT = 30000
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



    # Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")




    # Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space")
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")

    return beam, oe0


def run_shadow4_source(oe0):

    # beam.genSource(oe0)
    from shadow4.compatibility.translate_start_files import start00_to_bending_magnet
    from shadow4.compatibility.beam3 import Beam3
    import time

    bm4 = start00_to_bending_magnet(oe0)

    print(bm4.info())


    print(">>>>calling calculate_rays...")
    t00 = time.time()
    rays = bm4.calculate_rays(F_COHER=oe0.F_COHER,NRAYS=oe0.NPOINT,SEED=oe0.ISTAR1,
                                 EPSI_DX=oe0.EPSI_DX,EPSI_DZ=oe0.EPSI_DZ,verbose=False)
    t11 = time.time()-t00
    print(">>>> time for %d rays: %f s, %f min, "%(oe0.NPOINT,t11,t11/60))
    beam = Beam3.initialize_from_array(rays)

    return beam,oe0

if __name__ == "__main__":

    from srxraylib.plot.gol import plot, set_qt

    set_qt()

    beam_3, oe0_3 = run_shadow3_source()

    beam_4, oe0_3 = run_shadow4_source(oe0_3)


    Shadow.ShadowTools.plotxy(beam_3,1,4,nbins=101,nolost=1,title="Phase space X shadow3")
    Shadow.ShadowTools.plotxy(beam_4,1,4,nbins=101,nolost=1,title="Phase space X shadow4")

    Shadow.ShadowTools.plotxy(beam_3,3,6,nbins=101,nolost=1,title="Phase space Z shadow3")
    Shadow.ShadowTools.plotxy(beam_4,3,6,nbins=101,nolost=1,title="Phase space Z shadow4")

    Shadow.ShadowTools.histo1(beam_3,11,nbins=101,nolost=1)
    Shadow.ShadowTools.histo1(beam_4,11,nbins=101,nolost=1)

    print("Identical? ",beam_4.identical(beam_3))

    if beam_4.identical(beam_3):
        print(">>> BM sources with shadow4 and shadow3 are identical")
    else:
        beam_4.difference(beam_3)

