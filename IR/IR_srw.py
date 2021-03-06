try:
    from oasys_srw.srwlib import *
    from oasys_srw.uti_plot import *
except:
    from srwlib import *
    from uti_plot import *

import numpy

if not srwl_uti_proc_is_master(): exit()

# # //Electron Beam parameters
# ElecEnergy= 2.0    # // in GeV for ALR storage ring
# ElecCurrent = 0.50 # //[A] Electron Current
#
# EnergySpread=0.00095
# EmittanceX=2.0
# EmittanceZ=0.04
# BetaX=0.34
# BetaZ=24.26
# AlphaX=0.827
# AlphaZ=-10.7
# Dispers_X=0.031 #//(m)
# Dispers_Deriv=-0.06  #//(r)  ..
#
# # // Magnetic Field parameters
# FieldValue= 0.87            #// in Tesla for the magnt at ALS
# FieldCenter= 0              #// (m) Center of magnetic field
# FieldLength= 2.5            #//.5 //(m) length of the magnetic field
# FieldPointsNumber= 50000    #// Number of points for calculation ****************
# MagFieldLength=0.5          #//(m)
# MagSourcePoint=0            #//.269489//.225 //(m)
# MagFringeField=100          #//(mm)



##########################
distance = 2.935
max_div_h = 0.1
max_div_v = 0.015
##########################

try:
    from oasys_srw.srwlib import *
    from oasys_srw.uti_plot import *
except:
    from srwlib import *
    from uti_plot import *

import numpy

if not srwl_uti_proc_is_master(): exit()

####################################################
# LIGHT SOURCE

part_beam = SRWLPartBeam()
part_beam.Iavg               = 0.5
part_beam.partStatMom1.x     = 0.0
part_beam.partStatMom1.y     = 0.0
part_beam.partStatMom1.z     = 0.0
part_beam.partStatMom1.xp    = 0.0
part_beam.partStatMom1.yp    = 0.0
part_beam.partStatMom1.gamma = 3913.9027331571065
part_beam.arStatMom2[0]      = (7e-6)**2  #8.911024999999999e-10
part_beam.arStatMom2[1]      = (0e-6)**2  #-5.789e-11
part_beam.arStatMom2[2]      = (10e-6)**2  #3.466912647058823e-10
part_beam.arStatMom2[3]      = (10e-6)**2  #1.6982000000000002e-09
part_beam.arStatMom2[4]      = (0e-6)**2  #7.49e-10
part_beam.arStatMom2[5]      = (7e-6)**2  #3.582235779060181e-09
part_beam.arStatMom2[10]     = 9.025e-07


# SrwMagDipole("MagALSBZ_fld",1,-1.005, 0.5, 10, 0.876)     #        //Dipole 6
# SrwMagDipole("MagALSBZ_fld",2,-0.5275,0.305,10,-.16)      #  //RB6-7
# SrwMagDipole("MagALSBZ_fld",2,-.05,0.5, 10, 0.876)        #       ///Dipole 7
# SrwMagDipole("MagALSBZ_fld",2,0.4275, 0.305, 10, -0.16)   # //RB7-8
# SrwMagDipole("MagALSBZ_fld",2, 0.905,0.5, 10, 0.849)      #    //Dipole8

magnetic_structure1 = SRWLMagFldM(_G=-0.16, _m=1, _n_or_s='n', _Leff=0.305)
magnetic_structure2 = SRWLMagFldM(_G=0.867,    _m=1, _n_or_s='n', _Leff=0.500)
magnetic_structure3 = SRWLMagFldM(_G=-0.16, _m=1, _n_or_s='n', _Leff=0.305)
magnetic_structure4 = SRWLMagFldM(_G=0.849, _m=1, _n_or_s='n', _Leff=0.500)


magnetic_field_container = SRWLMagFldC(_arMagFld=[magnetic_structure1,magnetic_structure2,magnetic_structure3,magnetic_structure4],
                                    _arXc=array('d', [0.0, 0.0, 0.0, 0.0]),
                                    _arYc=array('d', [0.0, 0.0, 0.0, 0.0]),
                                    _arZc=array('d', [-0.5, 0.0, 0.5, 1.0]))

mesh = SRWLRadMesh(_eStart=0.1,
                   _eFin  =0.1,
                   _ne    =1,
                   _xStart= - max_div_h * distance, # -1.25,
                   _xFin  =   max_div_h * distance, # 1.25, #0.5,
                   _nx    =300,
                   _yStart= - max_div_v * distance, # -0.150,
                   _yFin  =   max_div_v * distance, # 0.150,
                   _ny    =200,
                   _zStart=distance) # 10.0)

stk = SRWLStokes()
stk.allocate(1,300,200)
stk.mesh = mesh

wfr = SRWLWfr()
wfr.allocate(mesh.ne, mesh.nx, mesh.ny)
wfr.mesh = mesh
wfr.partBeam = part_beam

initial_mesh = deepcopy(wfr.mesh)
# srwl.CalcElecFieldSR(wfr, 0, magnetic_field_container, [2,0.001,0.0,0.0,20000,1,0.0])
srwl.CalcElecFieldSR(wfr, 0, magnetic_field_container, [2,0.001,-0.75,1.25,20000,1,0.0])


print(">>>>>>>>>>>>>>>>>>>>>>>",wfr.Rx,wfr.Ry)



#
# native plots
#
if False:
    mesh0 = deepcopy(wfr.mesh)
    arI = array('f', [0]*mesh0.nx*mesh0.ny)
    srwl.CalcIntFromElecField(arI, wfr, 6, 0, 3, mesh0.eStart, 0, 0)
    arIx = array('f', [0]*mesh0.nx)
    srwl.CalcIntFromElecField(arIx, wfr, 6, 0, 1, mesh0.eStart, 0, 0)
    arIy = array('f', [0]*mesh0.ny)
    srwl.CalcIntFromElecField(arIy, wfr, 6, 0, 2, mesh0.eStart, 0, 0)
    #save ascii file with intensity
    #srwl_uti_save_intens_ascii(arI, mesh0, <file_path>)
    plotMesh0x = [1000*mesh0.xStart, 1000*mesh0.xFin, mesh0.nx]
    plotMesh0y = [1000*mesh0.yStart, 1000*mesh0.yFin, mesh0.ny]
    uti_plot2d1d (arI, plotMesh0x, plotMesh0y, labels=['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity Before Propagation'])
    uti_plot_show()

#
# oasys plots
#
from srxraylib.plot.gol import plot_image, plot
from wofrysrw.propagator.wavefront2D.srw_wavefront import SRWWavefront
import copy

self = wfr
w_srw_oasys = SRWWavefront(_arEx=copy.deepcopy(self.arEx),
                         _arEy=copy.deepcopy(self.arEy),
                         _typeE=self.numTypeElFld,
                         _eStart=self.mesh.eStart,
                         _eFin=self.mesh.eFin,
                         _ne=self.mesh.ne,
                         _xStart=self.mesh.xStart,
                         _xFin=self.mesh.xFin,
                         _nx=self.mesh.nx,
                         _yStart=self.mesh.yStart,
                         _yFin=self.mesh.yFin,
                         _ny=self.mesh.ny,
                         _zStart=self.mesh.zStart,
                         _partBeam=self.partBeam)



w_wofry = w_srw_oasys.toGenericWavefront()




x = w_wofry.get_coordinate_x()
y = w_wofry.get_coordinate_y()
intensity = w_wofry.get_intensity()
print("Total intensity: %g photons"%(intensity.sum()*1e6*(x[1]-x[0])*(y[1]-y[0])))


plot_image(intensity,x,y,show=0,aspect='auto')

plot(x/distance,intensity.sum(axis=1)*(y[1]*y[0]),xtitle="X' [rad]",ytitle = "photons/s/m",show=0)

plot(y/distance,intensity.sum(axis=0)*(x[1]*x[0]),xtitle="Y' [rad]",ytitle = "photons/s/m",show=1)


from wofrysrw.util.srw_hdf5 import save_wfr_2_hdf5
# save_wfr_2_hdf5(wfr,"/home/manuel/Oasys/wfr_ir_0p1_all_magnets.h5",subgroupname="wfr",intensity=True,phase=True,overwrite=True)
# save_wfr_2_hdf5(wfr,"/home/manuel/Oasys/wfr_ir_0p1_Mag7.h5",subgroupname="wfr",intensity=True,phase=True,overwrite=True)
