import Shadow
import numpy
def run_shadow(reflection=b'/home/manuel/Oasys/Ge220_3_25'):

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
    oe1 = Shadow.OE()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.FDISTR = 1
    oe0.FSOUR = 1
    oe0.F_COLOR = 3
    oe0.F_PHOT = 0
    oe0.F_POLAR = 0
    oe0.HDIV1 = 0.05
    oe0.HDIV2 = 0.05
    oe0.IDO_VX = 0
    oe0.IDO_VZ = 0
    oe0.IDO_X_S = 0
    oe0.IDO_Y_S = 0
    oe0.IDO_Z_S = 0
    oe0.ISTAR1 = 5676561
    oe0.NPOINT = 100000
    oe0.PH1 = 4995.0
    oe0.PH2 = 5005.0
    oe0.VDIV1 = 0.05
    oe0.VDIV2 = 0.05
    oe0.WXSOU = 2.0
    oe0.WZSOU = 0.02

    oe1.DUMMY = 0.1
    oe1.FCYL = 1
    oe1.FHIT_C = 1
    # oe1.FILE_REFL = b'/home/manuel/Oasys/Ge220_3_25'
    oe1.FILE_REFL = reflection
    oe1.FMIRR = 1
    oe1.FWRITE = 1
    oe1.F_CENTRAL = 1
    oe1.F_CRYSTAL = 1
    oe1.F_EXT = 1
    oe1.F_JOHANSSON = 1
    oe1.PHOT_CENT = 5000.0
    oe1.RLEN1 = 25.0
    oe1.RLEN2 = 25.0
    oe1.RMIRR = 250.0
    oe1.RWIDX1 = 7.0
    oe1.RWIDX2 = 7.0
    oe1.R_JOHANSSON = 500.0
    oe1.R_LAMBDA = 5000.0
    oe1.T_IMAGE = 309.917
    oe1.T_INCIDENCE = 51.6891770286
    oe1.T_REFLECTION = 51.6891770286
    oe1.T_SOURCE = 309.917



    #Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")


    #
    #run optical element 1
    #
    print("    Running optical element: %d"%(1))
    if iwrite:
        oe1.write("start.01")

    beam.traceOE(oe1,1)

    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")


    return beam


# reflection = b'/home/manuel/Oasys/Ge111_3_25'
reflection = b'/home/manuel/Oasys/Ge220_3_25'
beam = run_shadow(reflection=reflection)

# https://github.com/oasys-kit/shadow3/blob/master/Shadow/ShadowLibExtensions.py
print("intensity: ",beam.intensity(nolost=1))

tkt = beam.histo2(1,3,nolost=1,ref=23,nbins=101,calculate_widths=1)

for key in tkt.keys():
    print(key)

print(tkt["histogram"].shape,tkt["bin_h_center"].shape,tkt["bin_v_center"].shape)

from srxraylib.plot.gol import plot_image, plot

plot_image(tkt["histogram"],tkt["bin_h_center"],tkt["bin_v_center"],aspect='auto')

print(beam.rays.shape)


x = beam.getshonecol(1)
energy = beam.getshonecol(11)
print(energy)
# plot(x,energy)


Shadow.ShadowTools.plotxy(beam,1,3,nbins=101,nolost=1,title="Real space")
Shadow.ShadowTools.histo1(beam,11,nolost=1,ref=23)
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")