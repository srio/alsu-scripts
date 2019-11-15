

import scipy.constants as codata
import numpy
from srxraylib.plot.gol import plot, set_qt
from srxraylib.sources.srfunc import sync_ang

def eV_to_m(energy_eV):
    return codata.h * codata.c / codata.e / energy_eV

def eV_to_um(energy_eV):
    return 1e6*eV_to_m(energy_eV)

def eV_to_inversecm(energy_eV):
    wav = eV_to_m(energy_eV)
    return 1./(100*wav)

def flux_vs_psi(divergence = numpy.linspace(-50,50,1000),electronEnergy = 2.0,r_m = 7.6136,
                photon_energy=0.1,i_a=0.5):



    codata_mee = 1e-6 * codata.m_e * codata.c ** 2 / codata.e
    cte = (3.0e0/4/numpy.pi) * codata.h * codata.c * numpy.power(1e3 / codata_mee,3) / codata.e


    ec_ev = cte * numpy.power(electronEnergy,3) / r_m


    f_1 = sync_ang(1,divergence,polarization=0, \
        e_gev=electronEnergy,i_a=i_a,hdiv_mrad=1.0,r_m=r_m,energy=photon_energy,ec_ev=ec_ev)

    return f_1

if __name__ == "__main__":
    set_qt()

    print("12000 eV = %f A"%(1e10*eV_to_m(12000)))

    print("1 eV = %f um" % eV_to_um(1.))
    print(".1 eV = %f um" % eV_to_um(.1))
    print(".01 eV = %f um" % eV_to_um(.01))

    print("1 eV = %f cm-1" %   eV_to_inversecm(1.))
    print(".1 eV = %f cm-1" %  eV_to_inversecm(.1))
    print(".01 eV = %f cm-1" % eV_to_inversecm(.01))

    divergence = numpy.linspace(0,100,101)


    wavelength = eV_to_m(0.01)
    FWHM_p01 = 0.885 * wavelength / (divergence *1e-3)
    wavelength = eV_to_m(0.1)
    FWHM_p1 = 0.885 * wavelength / (divergence *1e-3)
    wavelength = eV_to_m(1)
    FWHM_1 = 0.885 * wavelength / (divergence *1e-3)

    plot(divergence,FWHM_1*1e3,
         divergence, FWHM_p1 * 1e3,
         divergence, FWHM_p01 * 1e3,
         ytitle="FWHM [mm]",xtitle="Divergence [mrad]",yrange=[0,10],
         legend=["E=1 eV","E=0.1 eV","E=0.01 eV"])

    divergence2 = numpy.linspace(-30,30,1000)

    f_1   = flux_vs_psi(divergence=divergence2,electronEnergy=2.0,r_m=7.6136,photon_energy=1)
    f_p1  = flux_vs_psi(divergence=divergence2,electronEnergy=2.0,r_m=7.6136,photon_energy=0.1)
    f_p01 = flux_vs_psi(divergence=divergence2,electronEnergy=2.0,r_m=7.6136,photon_energy=0.01)


    plot(divergence2, f_1/f_1.max(),
         divergence2, f_p1/f_p1.max(),
         divergence2, f_p01/f_p01.max(),
         ytitle="Normalized Flux",xtitle="Psi [mrad]",
         legend=["E=1 eV","E=0.1 eV","E=0.01 eV"])