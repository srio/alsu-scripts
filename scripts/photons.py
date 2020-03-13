# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 13:13:57 2019

@author: BRUMUND
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 09:25:10 2019

@author: brumund
"""

#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as codata
m2ev = codata.c * codata.h / codata.e      # lambda(m)  = m2eV / energy(eV)


# def run():
#   print("Starting calculation for photon distribution")
#   E=np.arange(12000,20001,1000)
#   fwhm=np.zeros((E.size,2))
#
#   for Ei in E:
#     #print("" + str())
#     progr=round(E.searchsorted(Ei)/E.size*100)
#     print("Progress: " + str(progr)+ "%" )
#     beam,_,fwhm_h,fwhm_v=ray_tracing(Ei,0.1)
#     print(">>>>>>>",Ei,fwhm_h,fwhm_v)
#     fwhm[E.searchsorted(Ei)]=[fwhm_h,fwhm_v]
#
#   return E,fwhm,beam
    



def get_sigmas_radiation(photon_energy,undulator_length):
    lambdan = m2ev / photon_energy
    return 2.740/4/np.pi*np.sqrt(lambdan*undulator_length),0.69*np.sqrt(lambdan/undulator_length),

def get_sigmas_EBS(photon_energy,undulator_length):
    sr, srp = get_sigmas_radiation(photon_energy, undulator_length)

    sx, sz, sxp, szp = 3.01836e-05, 3.63641e-06, 4.36821e-06, 1.37498e-06

    Sx = np.sqrt(sx ** 2 + sr ** 2)
    Sz = np.sqrt(sz ** 2 + sr ** 2)
    Sxp = np.sqrt(sxp ** 2 + srp ** 2)
    Szp = np.sqrt(szp ** 2 + srp ** 2)

    return Sx,Sz,Sxp,Szp


# Script output from Shadow Info widget(python script)
def ray_tracing(E,BW):
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
    
    #
    # Define variables. See meaning of variables in: 
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml 
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    NEWSIGMAS = get_sigmas_EBS(E,0.016*125)

    oe0.FDISTR = 3
    oe0.FSOUR = 3
    oe0.F_COLOR = 3
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.0
    oe0.HDIV2 = 0.0
    oe0.ISTAR1 = 5676561
    oe0.NPOINT = 50000
    oe0.PH1 = E
    oe0.PH2 = E
    oe0.SIGDIX = NEWSIGMAS[0]
    oe0.SIGDIZ = NEWSIGMAS[1]
    oe0.SIGMAX = NEWSIGMAS[2]
    oe0.SIGMAZ = NEWSIGMAS[3]
    oe0.VDIV1 = 0.0
    oe0.VDIV2 = 0.0
    
    oe1.DUMMY = 100.0
    oe1.FWRITE = 0
    oe1.F_REFRAC = 2
    oe1.F_SCREEN = 1
    oe1.N_SCREEN = 1
    oe1.T_IMAGE = 0.0
    oe1.T_INCIDENCE = 0.0
    oe1.T_REFLECTION = 180.0
    oe1.T_SOURCE = 30.0
    
    
    
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
    
    
    # Shadow.ShadowTools.plotxy(beam,1,3,nbins=101,nolost=1,title="Real space E=%f"%E)
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")
    
    
    hist2=beam.histo2(1,3,nbins=50,nbins_h=100,nbins_v=100)
    fwhm_h=hist2.get("fwhm_h")
    fwhm_v=hist2.get("fwhm_v")
    
    return beam,E,fwhm_h,fwhm_v
    
    #Shadow.ShadowTools.plotxy(beam,1,3,nbins=101,nolost=1,title="Real space")
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")
    
  
if __name__== "__main__":

  print("Starting calculation for photon distribution")
  E = np.linspace(12000,50000,5) # np.arange(12000, 20001, 1000)
  fwhm = np.zeros((E.size, 2))

  for Ei in E:
      progr = round(E.searchsorted(Ei) / E.size * 100)
      print("Progress: " + str(progr) + "%")
      beam, _, fwhm_h, fwhm_v = ray_tracing(Ei, 0.1)
      print(">>>>>>>", Ei, fwhm_h, fwhm_v)
      fwhm[E.searchsorted(Ei)] = [fwhm_h, fwhm_v]


  plt.plot(E,fwhm[:,0],)
  plt.title("Horizontal")
  plt.show()
  plt.plot(E, fwhm[:, 1])
  plt.title("Vertical")
  plt.show()
  
  
#  beam,E,fwhm_h,fwhm_v=ray_tracing(30000,1)
#  Shadow.ShadowTools.plotxy(beam,1,3,nbins=101,nolost=1,title="Real space")