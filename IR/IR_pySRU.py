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


    plot(z*codata.c,x0,z*codata.c,x,legend=["sinusoidal","modified"])

    simulation_test_input.trajectory = Trajectory(t=t,x=x,y=y,z=z,v_x=v_x,v_y=v_y,v_z=v_z,
                                  a_x=a_x,a_y=a_y,a_z=a_z)



    return simulation_test



if __name__ == "__main__":
    from srxraylib.plot.gol import set_qt

    set_qt()


    print("======================================================================")
    print("======      Undulator U18 from ESRF with K=1.68                =======")
    print("======================================================================")

    beam_ALSU = ElectronBeam(Electron_energy=2.0, I_current=0.5)
    # ESRF18 = Undulator(K=1.68, period_length=0.018, length=2.0)
    id = Undulator(K=2.0, period_length=0.4, length=1.6)

    #
    # radiation in a defined mesh with only one point to save time
    #

    distance = None
    simulation_test = create_simulation(magnetic_structure=id,
                                          electron_beam=beam_ALSU,
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
    #     shift = 0.0
    #     omega0 = simulation_test.radiation_fact.photon_frequency
    #     omega = omega0 + omega0*shift
    #     simulation_test.radiation_fact.photon_frequency = omega
    #
    #     print(">>>>>>>>>> E shifted: ",simulation_test.radiation_fact.energy_eV())
        simulation_test.change_XY_radiation(X=X, Y=Y, update_radiation=True)
    #
        plot(
            1e6 * simulation_test0.radiation.Y[0, :],
            simulation_test0.radiation.intensity[0, :],
            1e6*simulation_test.radiation.Y[0,:],
            simulation_test.radiation.intensity[0,:],
            xtitle="Y [urad]",ytitle="Intensity",
            legend=["sinusoidal","modified"],
            show=True)
    #
    # simulation_test.print_parameters()
    # # simulation_test.radiation.plot(
    # #     title=("radiation in a screen including central cone, harmonic %d,at D=" + repr(distance)) % harmonic_number)
    #
    # print(simulation_test.radiation.intensity.shape,simulation_test.radiation.X.shape,simulation_test.radiation.Y.shape)

