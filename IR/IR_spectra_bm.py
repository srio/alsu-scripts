import numpy
from srxraylib.plot.gol import plot, set_qt
from scipy.ndimage import gaussian_filter1d
import scipy.constants as codata
import srxraylib.sources.srfunc as srfunc
from IR_magnetic_field import get_magnetic_field_ALS,get_magnetic_field_ALSU



if __name__ == "__main__":

    from orangecontrib.xoppy.util.xoppy_bm_wiggler import xoppy_calc_bm
    set_qt()

    #
    # spectrum
    #

    # TYPE_CALC:
    # 0: 'Energy or Power spectra'
    # 1: 'Angular distribution (all wavelengths)'
    # 2: 'Angular distribution (one wavelength)'
    # 3: '2D flux and power (angular,energy) distribution'
    #
    # Radius
    # M1: 7.615618611829955
    # Radius
    # AB: 41.695511899769
    # Radius
    # M2: 7.857811429874018
    # Half - Divergence
    # M1: 0.032827274151
    # Half - Divergence
    # AB: 0.0038973019540000007
    # Half - Divergence
    # M2: 0.031841706445325


    npoint = 500
    e0 = numpy.zeros(500)
    out = numpy.zeros((npoint,4))

    MACHINE_NAME    = ["ALS","ALSU-M7","ALSU-AB","ALSU-M8"]
    BFIELD          = [1.27,0.87,0.16,0.87]
    BEAM_ENERGY_GEV = [1.9,2.0,2.0,2.0]



    #
    # full emission
    #


    HOR_DIV_MRAD = [69, 2e3 * 0.0328, 2e3 * 0.003897]  # not used if VER_DIV = 0
    # HOR_DIV_MRAD    = [1.0]*3

    VER_DIV_MRAD    = [1.0, 1.0, 1.0]


    for i in range(3):
        a6_T, fm, a, energy_ev = xoppy_calc_bm(
            TYPE_CALC=0,
            MACHINE_NAME=MACHINE_NAME[i],
            RB_CHOICE=1,
            MACHINE_R_M=0.0,
            BFIELD_T=BFIELD[i],
            BEAM_ENERGY_GEV=BEAM_ENERGY_GEV[i],
            CURRENT_A=0.5,
            HOR_DIV_MRAD=HOR_DIV_MRAD[i],
            VER_DIV=0, #2, #0,            <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            PHOT_ENERGY_MIN=0.001,
            PHOT_ENERGY_MAX=10000.0,
            NPOINTS=500,
            LOG_CHOICE=0,
            PSI_MRAD_PLOT=20.0,
            PSI_MIN = -VER_DIV_MRAD[i] * 0.5,
            PSI_MAX =  VER_DIV_MRAD[i] * 0.5,
            PSI_NPOINTS=500,
            FILE_DUMP=1)  # writes output to bm.spec

        e0 = a6_T[:, 0].copy()
        out[:,i] = a6_T[:, 5].copy()

    # print(a6_T.shape)



    # plot(e0, 2*out[:, 1] + out[:, 2],
    #      e0, out[:, 1],
    #      e0, out[:, 2],
    #      e0, out[:, 0],
    #      xlog=True, ylog=True, show=True, yrange=[1e11,1e16], #yrange=[1e5,1e10], #
    #      xtitle="Photon energy [eV]", ytitle="Flux [Photons/s/0.1%bw]", title="Flux (full emission)",
    #      legend=["ALSU: Mag7+Antibend+Mag8","ALSU: Mag7","ALSU: Antibend","ALS"],
    #      linestyle=["solid","solid","solid","dashed"])



    #
    # emission in accepted window
    #


    HOR_DIV_MRAD = [60, 66, 7.5, 29.5]  # not used if VER_DIV = 0
    VER_DIV_MRAD = [17, 14, 28,  34.0]


    for i in range(len(MACHINE_NAME)):
        a6_T, fm, a, energy_ev = xoppy_calc_bm(
            TYPE_CALC=0,
            MACHINE_NAME=MACHINE_NAME[i],
            RB_CHOICE=1,
            MACHINE_R_M=0.0,
            BFIELD_T=BFIELD[i],
            BEAM_ENERGY_GEV=BEAM_ENERGY_GEV[i],
            CURRENT_A=0.5,
            HOR_DIV_MRAD=HOR_DIV_MRAD[i],
            VER_DIV=2, #0,            <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            PHOT_ENERGY_MIN=0.001,
            PHOT_ENERGY_MAX=10000.0,
            NPOINTS=500,
            LOG_CHOICE=0,
            PSI_MRAD_PLOT=20.0,
            PSI_MIN = -VER_DIV_MRAD[i] * 0.5,
            PSI_MAX =  VER_DIV_MRAD[i] * 0.5,
            PSI_NPOINTS=500,
            FILE_DUMP=1)  # writes output to bm.spec

        e0 = a6_T[:, 0].copy()
        out[:,i] = a6_T[:, 5].copy()




    plot(e0, out[:, 1] + out[:, 2] + out[:, 3],
         e0, out[:, 1],
         e0, out[:, 2],
         e0, out[:, 0],
         e0, out[:, 3],
         xlog=True, ylog=True, show=True, xrange=[1e-3,1e1], yrange=[1e12,1e15], #yrange=[1e5,1e10], #
         xtitle="Photon energy [eV]", ytitle="Flux [Photons/s/0.1%bw]", title="Flux (room constrainted emission)",
         legend=["ALSU: Mag7+Antibend+Mag8","ALSU: Mag7","ALSU: Antibend","ALS", "ALSU: Mag8"],
         linestyle=["solid","solid","solid","dashed","solid"])
