from pySRU.ElectronBeam import ElectronBeam
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
from pySRU.Simulation import create_simulation, Simulation
from pySRU.TrajectoryFactory import TRAJECTORY_METHOD_ANALYTIC,TRAJECTORY_METHOD_ODE
from pySRU.RadiationFactory import RADIATION_METHOD_NEAR_FIELD, RADIATION_METHOD_APPROX_FARFIELD

import numpy as np



from pySRU.ElectronBeam import ElectronBeam
import numpy as np
import scipy.stats
from  scipy.constants import physical_constants
#from mpl_toolkits.mplot3d import Axes3D
import scipy.constants as codata
import time
from scipy.interpolate import interp1d
from pySRU.Radiation import Radiation
from pySRU.MagneticField import MagneticField
from pySRU.SourceUndulatorPlane import SourceUndulatorPlane
from pySRU.SourceBendingmagnet import SourceBendingMagnet
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
from pySRU.MagneticStructureBendingMagnet import MagneticStructureBendingMagnet as BM

from pySRU.Trajectory import Trajectory

from pySRU.TrajectoryFactory import TrajectoryFactory, TRAJECTORY_METHOD_ANALYTIC,TRAJECTORY_METHOD_ODE

from pySRU.RadiationFactory import RadiationFactory,RADIATION_METHOD_NEAR_FIELD, \
                                 RADIATION_METHOD_APPROX_FARFIELD


eV_to_J = codata.e

from srxraylib.plot.gol import plot_image, plot

#
# def create_simulation11(magnetic_structure,electron_beam, magnetic_field=None, photon_energy=None,
#                       traj_method=TRAJECTORY_METHOD_ANALYTIC,Nb_pts_trajectory=None,
#                       rad_method=RADIATION_METHOD_APPROX_FARFIELD, Nb_pts_radiation=101,
#                       initial_condition=None, distance=None,XY_are_list=False,X=None,Y=None) :
#
#     if type(magnetic_structure)==Undulator :
#         source=SourceUndulatorPlane(undulator=magnetic_structure,
#                                     electron_beam=electron_beam, magnetic_field=magnetic_field)
#         print("Calculating undulator source...")
#     elif type(magnetic_structure)==BM:
#         source = SourceBendingMagnet(magnetic_structure=magnetic_structure,
#                                       electron_beam=electron_beam, magnetic_field=magnetic_field)
#         print("Calculating bending magnet source...")
#     else :
#         raise Exception('magnet type unknown')
#
#     if photon_energy==None :
#         omega=source.choose_photon_frequency()
#     else :
#         omega = photon_energy * eV_to_J / codata.hbar
#
#     #
#     # XY grid
#     #
#     if distance==None and (rad_method==RADIATION_METHOD_NEAR_FIELD) :
#         distance=source.choose_distance_automatic(2)
#
#     if Nb_pts_trajectory==None :
#         Nb_pts_trajectory = int(source.choose_nb_pts_trajectory(0.01,photon_frequency=omega))
#
#     if X is None or Y is None :
#         if (X != None) :
#             Y=X
#         elif Y != None :
#             X=Y
#         else :
#             theta_max=source.choose_angle_deflection_max()
#             if distance==None :
#                 X_max=theta_max
#                 Y_max=theta_max
#             else :
#                 X_max = distance * theta_max
#                 Y_max = distance * theta_max
#             X = np.linspace(0.0, X_max, Nb_pts_radiation)
#             Y = np.linspace(0.0, Y_max, Nb_pts_radiation)
#
#     if type(X) == float:
#         X= np.linspace(0.0, X, Nb_pts_radiation)
#     if type(Y) == float:
#         Y = np.linspace(0.0, Y, Nb_pts_radiation)
#
#     # if X.shape != Y.shape :
#     #     raise Exception('X and Y must have the same shape')
#     Nb_pts_radiation = X.size # len(X.flatten())
#
#     #print('step 1')
#     traj_fact=TrajectoryFactory(Nb_pts=Nb_pts_trajectory,method=traj_method,initial_condition=initial_condition)
#
#     if (traj_fact.initial_condition == None):
#         # print('crearte initial cond automat')
#         traj_fact.initial_condition = source.choose_initial_contidion_automatic()
#
#
#
#     #print('step 2')
#
#     rad_fact=RadiationFactory(method=rad_method, photon_frequency=omega)
#
#
#     #print('step 3')
#     trajectory=traj_fact.create_from_source(source=source)
#     #print('step 4')
#     radiation = rad_fact.create_for_one_relativistic_electron(trajectory=trajectory, source=source, XY_are_list=XY_are_list,
#                                                               distance=distance, X=X, Y=Y)
#
#     #print('step 5')
#     return Simulation(source=source, trajectory_fact=traj_fact,
#                                radiation_fact=rad_fact, trajectory=trajectory, radiation=radiation)
#
#
#
#
#

def change_trajectory(simulation_test_input):
    simulation_test = simulation_test_input.copy()

    t = simulation_test.trajectory.t

    x = simulation_test.trajectory.x
    y = simulation_test.trajectory.y
    z = simulation_test.trajectory.z


    x0 = x.copy()

    x *= codata.c
    xM = x.max()
    x /= xM


    # x =  np.sign(x)
    a = 0.099/2
    tt = t*codata.c + 20
    x = (2/a) * (tt - a * np.fix(tt/a + 0.5) ) * (-1)**(np.fix(tt/a-0.5))

    x = x * xM / codata.c






    # v_x = simulation_test.trajectory.v_x * 0
    # v_y = simulation_test.trajectory.v_y * 0
    # v_z = simulation_test.trajectory.v_z * 0


    v_x = np.gradient(x, t)
    v_y = np.gradient(y, t)
    v_z = np.gradient(z, t)


    a_x = simulation_test.trajectory.a_x * 0
    a_y = simulation_test.trajectory.a_y * 0
    a_z = simulation_test.trajectory.a_z * 0


    plot(z,x0,z,x,legend=["sinusoidal","modified"])

    simulation_test_input.trajectory = Trajectory(t=t,x=x,y=y,z=z,v_x=v_x,v_y=v_y,v_z=v_z,
                                  a_x=a_x,a_y=a_y,a_z=a_z)



    return simulation_test


def gaussian_fit(y,x):

    #
    # Fitting
    #
    from silx.math.fit.functions import sum_gauss
    from silx.math.fit import fittheories
    from silx.math.fit.fitmanager import FitManager

    p = [y.max(),0,5]

    fit = FitManager()
    fit.setdata(x=x, y=y)
    fit.loadtheories(fittheories)
    fit.settheory('Gaussians')
    fit.estimate()
    fit.runfit()


    print("Searched parameters = %s" % p)
    print("Obtained parameters : ")
    dummy_list = []
    for param in fit.fit_results:
        print(param['name'], ' = ', param['fitresult'])
        dummy_list.append(param['fitresult'])
    print("chisq = ", fit.chisq)
    fwhm_txt = "FWHM of fit = %5.3f um"%(fit.fit_results[2]['fitresult'])
    y1 = sum_gauss(x, *dummy_list)

    return y1

def gaussian_fit_scan(simulation_test0,shift_min=-5e-2,shift_max=5e-2,shift_n=5,do_plot=True):
    from silx.math.fit.functions import sum_gauss
    from silx.math.fit import fittheories
    from silx.math.fit.fitmanager import FitManager

    SHIFT = np.linspace(shift_min,shift_max,shift_n)





    CHI = np.zeros_like(SHIFT)
    KUR = np.zeros_like(SHIFT)

    for i,shift in enumerate(SHIFT):

        simulation_test = simulation_test0.copy()
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>",SHIFT[i])

        omega0 = simulation_test.radiation_fact.photon_frequency
        omega = omega0 + omega0 * shift
        simulation_test.radiation_fact.photon_frequency = omega

        simulation_test.update_radiation()

        #
        #
        #

        x = 1e6 * simulation_test.radiation.Y[0, :]
        y = simulation_test.radiation.intensity[0, :]

        y /= y.max()

        #
        # Fitting
        #

        p = [y.max(),0,0.3*x.max()]

        fit = FitManager()
        fit.setdata(x=x, y=y)
        fit.loadtheories(fittheories)
        fit.settheory('Gaussians')
        fit.estimate()
        fit.runfit()


        print("Searched parameters = %s" % p)
        print("Obtained parameters : ")
        dummy_list = []
        for param in fit.fit_results:
            print(param['name'], ' = ', param['fitresult'])
            dummy_list.append(param['fitresult'])
        print("chisq = ", fit.chisq)
        fwhm_txt = "FWHM of fit = %5.3f um"%(fit.fit_results[2]['fitresult'])
        y1 = sum_gauss(x, *dummy_list)

        # plot(x,y,x,y1,legend=["data","fit"])

        print(">>>> Kurtosis0: ", scipy.stats.kurtosis(y,axis=None))
        print(">>>> Kurtosis: ", scipy.stats.kurtosis(y1,axis=None))

        print(">>>>>>> chi**2", ((y-y1)**2).sum())
        KUR[i] =  scipy.stats.kurtosis(y1,axis=None)
        CHI[i] = ((y-y1)**2).sum()

    plot(SHIFT,CHI,title="CHI",show=False)
    plot(SHIFT,KUR,title="kur")




if __name__ == "__main__":
    from srxraylib.plot.gol import set_qt

    set_qt()

    #
    # print("======================================================================")
    # print("======      Undulator from X-ray data booklet                  =======")
    # print("====== fig 2.5 in  http://xdb.lbl.gov/Section2/Sec_2-1.html    =======")
    # print("======================================================================")
    #
    # # note that the flux in the reference fig 2.6 is a factor 10 smaller than the calculated here.
    # # This factor comes from the units:
    # #     here: phot / s  / A / 0.1%bw / (mrad)^2
    # #     ref : phot / s  / A /   1%bw / (0.1 mrad)^2
    #
    # undulator_test = Undulator(K=1.87, period_length=0.035, length=0.035 * 14)
    # electron_beam_test = ElectronBeam(Electron_energy=1.3, I_current=1.0)
    #
    # simulation_test = create_simulation(magnetic_structure=undulator_test,electron_beam=electron_beam_test,
    #                     magnetic_field=None, photon_energy=None,
    #                     traj_method=TRAJECTORY_METHOD_ANALYTIC,Nb_pts_trajectory=None,
    #                     rad_method=RADIATION_METHOD_APPROX_FARFIELD, Nb_pts_radiation=101,
    #                     initial_condition=None, distance=None,XY_are_list=False,X=None,Y=None)
    #
    #
    # simulation_test.print_parameters()
    #
    # simulation_test.trajectory.plot_3D(title="Electron Trajectory")
    #
    # simulation_test.radiation.plot(title="Flux in far field vs angle")


    print("======================================================================")
    print("======      Undulator U18 from ESRF with K=1.68                =======")
    print("======================================================================")

    beam_ESRF = ElectronBeam(Electron_energy=6.0, I_current=0.2)
    # ESRF18 = Undulator(K=1.68, period_length=0.018, length=2.0)
    ESRF18 = Undulator(K=1.0, period_length=0.1, length=1.0)

    #
    # radiation in a defined mesh with only one point to save time
    #

    distance = None
    simulation_test = create_simulation(magnetic_structure=ESRF18,
                                          electron_beam=beam_ESRF,
                                          magnetic_field=None,
                                          photon_energy=None,
                                          traj_method=TRAJECTORY_METHOD_ANALYTIC,
                                          Nb_pts_trajectory=None,
                                          rad_method=RADIATION_METHOD_APPROX_FARFIELD, Nb_pts_radiation=101,
                                          initial_condition=None,
                                          distance=distance,
                                          XY_are_list=False,
                                          X=np.array([0]),
                                          Y=np.array([0]))

    simulation_test0 = simulation_test.copy()
    #
    # central cone
    #

    if False:
        harmonic_number = 1
        print('create radiation in a screen including central cone for harmonic %d' % harmonic_number)
        simulation_test.calculate_on_central_cone(harmonic_number=harmonic_number, npoints_x=60, npoints_y=50)

        plot_image(simulation_test.radiation.intensity,
                   1e6 * simulation_test.radiation.X[:, 0],
                   1e6 * simulation_test.radiation.Y[0, :],
                   xtitle="X [urad]", ytitle="Y [urad]",
                   show=False)

    else:

        X = np.array([0])  # np.linspace(0.0, 12e-6, 50)
        Y = np.linspace(-75e-6, 75e-6, 500)
        # Y,X=np.meshgrid(X,Y)

        X_ones = np.ones_like(X)
        Y_ones = np.ones_like(Y)
        X = np.outer(X, Y_ones)
        Y = np.outer(X_ones, Y)

        change_trajectory(simulation_test)

        # simulation_test.trajectory.plot_2D()



        simulation_test0.change_XY_radiation(X=X, Y=Y, update_radiation=True)

        print(">>>>>>>>> E0: ", simulation_test.radiation_fact.energy_eV())
        shift = 0.0
        omega0 = simulation_test.radiation_fact.photon_frequency
        omega = omega0 + omega0*shift
        simulation_test.radiation_fact.photon_frequency = omega

        print(">>>>>>>>>> E shifted: ",simulation_test.radiation_fact.energy_eV())
        simulation_test.change_XY_radiation(X=X, Y=Y, update_radiation=True)

        plot(
            1e6 * simulation_test0.radiation.Y[0, :],
            simulation_test0.radiation.intensity[0, :],
            1e6*simulation_test.radiation.Y[0,:],
            simulation_test.radiation.intensity[0,:],
            xtitle="Y [urad]",ytitle="Intensity",
            legend=["sinusoidal","modified"],
            show=True)

    simulation_test.print_parameters()
    # simulation_test.radiation.plot(
    #     title=("radiation in a screen including central cone, harmonic %d,at D=" + repr(distance)) % harmonic_number)

    print(simulation_test.radiation.intensity.shape,simulation_test.radiation.X.shape,simulation_test.radiation.Y.shape)

    # simulation_test.trajectory.plot_2D()



    # y1 = gaussian_fit(1e6 * simulation_test.radiation.Y[0, :], 1e6 * simulation_test0.radiation.Y[0, :],)

    # plot(
    #     1e6 * simulation_test0.radiation.Y[0, :],
    #     simulation_test0.radiation.intensity[0, :],
    #     1e6 * simulation_test.radiation.Y[0, :],
    #     simulation_test.radiation.intensity[0, :],
    #     1e6 * simulation_test0.radiation.Y[0, :],
    #     y1,
    #     xtitle="Y [urad]", ytitle="Intensity",
    #     legend=["sinusoidal", "modified","fit sinusoidal"],
    #     show=True)

    gaussian_fit_scan(simulation_test0,shift_min=-10e-2,shift_max=10e-2,shift_n=25,do_plot=True)

    # plot(        1e6 * simulation_test0.radiation.Y[0, :],
    #     y1,)
    # #
    # # up to a given ring
    # #
    # harmonic_number = 1
    # ring_number = 1
    # print('create radiation a given screen including ring %d for harmonic %d' % (ring_number, harmonic_number))
    # simulation_test.calculate_until_ring_number(harmonic_number=harmonic_number, ring_number=ring_number, XY_are_list=False,
    #                                             npoints=51)
    # simulation_test.print_parameters()
    # simulation_test.radiation.plot(title=" radiation in a screen containing harmonic %d up to ring %d ring" %
#                                      (harmonic_number, ring_number))