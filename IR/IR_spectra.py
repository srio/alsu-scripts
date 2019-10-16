import numpy
from srxraylib.plot.gol import plot, set_qt
from scipy.ndimage import gaussian_filter1d
import scipy.constants as codata
import srxraylib.sources.srfunc as srfunc
from IR_magnetic_field import get_magnetic_field_ALS,get_magnetic_field_ALSU

def calculate_flux(y,B,M=0):
    # analyse M1

    B3 = B.copy()

    if M==1:
        ibad = numpy.argwhere(y > -0.2)
        B3[ibad] = 0.0
    elif M==2:
        ibad = numpy.argwhere( numpy.abs(y) > 0.2)
        B3[ibad] = 0.0
    elif M==3:
        ibad = numpy.argwhere(y < 0.2)
        B3[ibad] = 0.0
    else:
        pass


    tmp = numpy.vstack((y,B3)).T
    print(">>>>",tmp.shape,tmp[:,0])

    (traj, pars) = srfunc.wiggler_trajectory(
        b_from=1,
        inData=tmp, #"BM_first.b",
        nPer=1,
        nTrajPoints=501,
        ener_gev=1.9,
        per=0.01,
        kValue=1.0,
        trajFile="tmp.traj",
        shift_x_flag=5,
        shift_x_value=0.042,
        shift_betax_flag=4,
        shift_betax_value=0.0324)

    calculate_spectrum = True

    if calculate_spectrum:
        e, f, w = srfunc.wiggler_spectrum(traj,
                                          enerMin=0.0010,
                                          enerMax=10000.1,
                                          nPoints=500,
                                          electronCurrent=400.0 * 1e-3,
                                          outFile="spectrum.dat",
                                          elliptical=False)
        from srxraylib.plot.gol import plot

    return e,f,w





if __name__ == "__main__":
    set_qt()

    electron_energy_in_GeV = 1.9

    print("ALSU Radius M1: ", 1e9 / codata.c * electron_energy_in_GeV/0.876)
    print("ALSU Radius AB: ", 1e9 / codata.c * electron_energy_in_GeV/0.16)
    print("ALSU Radius M2: ", 1e9 / codata.c * electron_energy_in_GeV/0.8497)

    print("ALSU Half-Divergence M1: ", 0.5 * (0.500) / (1e9 / codata.c * electron_energy_in_GeV/0.876) )
    print("ALSU Half-Divergence AB: ", 0.5 * (0.325) / (1e9 / codata.c * electron_energy_in_GeV/0.16) )
    print("ALSU Half-Divergence M2: ", 0.5 * (0.500) / (1e9 / codata.c * electron_energy_in_GeV/0.8497) )



    #
    # spectrum
    #

    yy, B2 = get_magnetic_field_ALSU(do_plot=True)

    e1, f1, w1 = calculate_flux(yy, B2, M=1)
    e2, f2, w2 = calculate_flux(yy, B2, M=2)
    e3, f3, w3 = calculate_flux(yy, B2, M=3)
    e0, f0, w0 = calculate_flux(yy, B2, M=0)

    yy, B2 = get_magnetic_field_ALS(do_plot=True)
    eold, fold, wold = calculate_flux(yy, B2, M=0)

    plot(e0, f0,
         e1, f1,
         e2, f2,
         eold, fold,
         xlog=True, ylog=True, show=False, yrange=[1e11,1e16],
         xtitle="Photon energy [eV]", ytitle="Flux [Photons/s/0.1%bw]", title="Flux",
         legend=["ALSU: Mag7+Antibend+Mag8","ALSU: Mag7","ALSU: Antibend","ALS"],
         linestyle=["solid","solid","solid","dashed"])


    plot(e0, f0,
         e1, f1,
         e2, f2,
         eold, fold,
         xlog=True, ylog=False, show=True,yrange=[0,1e15],xrange=[1e-3,1],
         xtitle="Photon energy [eV]", ytitle="Flux [Photons/s/0.1%bw]", title="Flux",
         legend=["ALSU: Mag7+Antibend+Mag8","ALSU: Mag7","ALSU: Antibend","ALS"],
         linestyle=["solid","solid","solid","dashed"])