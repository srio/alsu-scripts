import numpy


def calculate_wavefront1D(wavelength=1e-10,
                          undulator_length=1.0, undulator_distance=10.0,
                          x_min=-0.1, x_max=0.1, number_of_points=101,
                          wavefront_position=0, add_random_phase=0):
    from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D

    sigma_r = 2.740 / 4 / numpy.pi * numpy.sqrt(wavelength * undulator_length)
    sigma_r_prime = 0.69 * numpy.sqrt(wavelength / undulator_length)

    wavefront1D = GenericWavefront1D.initialize_wavefront_from_range(x_min=x_min, x_max=x_max,
                                                                     number_of_points=number_of_points)
    wavefront1D.set_wavelength(wavelength)

    if wavefront_position == 0:  # Gaussian source
        wavefront1D.set_gaussian(sigma_x=sigma_r, amplitude=1.0, shift=0.0)
    elif wavefront_position == 1:  # Spherical source, Gaussian intensity
        wavefront1D.set_spherical_wave(radius=undulator_distance, center=0.0, complex_amplitude=complex(1, 0))
        # weight with Gaussian
        X = wavefront1D.get_abscissas()
        A = wavefront1D.get_complex_amplitude()
        sigma = undulator_distance * sigma_r_prime
        sigma_amplitude = sigma * numpy.sqrt(2)
        Gx = numpy.exp(-X * X / 2 / sigma_amplitude ** 2)
        wavefront1D.set_complex_amplitude(A * Gx)

    if add_random_phase:
        wavefront1D.add_phase_shifts(2 * numpy.pi * numpy.random.random(wavefront1D.size()))

    return wavefront1D


import numpy

from numba import jit, prange
@jit(nopython=True, parallel=True)
def goFromToSequential(field1, x1, y1, x2, y2, wavelength=1e-10, normalize_intensities=False):
    field2 = x2 * 0j # numpy.zeros_like(x2, dtype=complex)
    wavenumber = numpy.pi * 2 / wavelength

    for i in prange(field2.size):
        r = numpy.sqrt(numpy.power(x1 - x2[i], 2) + numpy.power(y1 - y2[i], 2))
        field2[i] = (field1 * numpy.exp(1.j * wavenumber * r)).sum()


    if normalize_intensities:
        field2 *= numpy.sqrt((numpy.abs(field1) ** 2).sum() / (numpy.abs(field2) ** 2).sum())
    return field2


def propagator1D_offaxis(input_wavefront, x2_oe, y2_oe, p, q, theta_grazing_in, theta_grazing_out=None,
                         zoom_factor=1.0, normalize_intensities=False):
    from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D

    if theta_grazing_out is None:
        theta_grazing_out = theta_grazing_in

    x1 = input_wavefront.get_abscissas()
    field1 = input_wavefront.get_complex_amplitude()
    wavelength = input_wavefront.get_wavelength()

    x1_oe = -p * numpy.cos(theta_grazing_in) + x1 * numpy.sin(theta_grazing_in)
    y1_oe = p * numpy.sin(theta_grazing_in) + x1 * numpy.cos(theta_grazing_in)

    # field2 is the electric field in the mirror
    field2 = goFromToSequential(field1, x1_oe, y1_oe, x2_oe, y2_oe,
                                wavelength=wavelength, normalize_intensities=normalize_intensities)

    x3 = x1 * zoom_factor

    x3_oe = q * numpy.cos(theta_grazing_out) + x3 * numpy.sin(theta_grazing_out)
    y3_oe = q * numpy.sin(theta_grazing_out) + x3 * numpy.cos(theta_grazing_out)

    # field2 is the electric field in the image plane
    field3 = goFromToSequential(field2, x2_oe, y2_oe, x3_oe, y3_oe,
                                wavelength=wavelength, normalize_intensities=normalize_intensities)

    output_wavefront = GenericWavefront1D.initialize_wavefront_from_arrays(x3, field3, wavelength=wavelength)

    return output_wavefront


def calculate_output_wavefront_after_grazing_reflector1D(input_wavefront, shape=1, radius=10000.0,
                                                         grazing_angle_in=1.5e-3, grazing_angle_out=1.5e-3,
                                                         p_distance=1.0,
                                                         q_distance=1.0,
                                                         zoom_factor=1.0,
                                                         error_flag=0, error_file="", write_profile=0):
    x1 = input_wavefront.get_abscissas()
    field1 = input_wavefront.get_complex_amplitude()

    if error_flag == 0:  # no profile file
        x2_oe = x1 / numpy.sin(grazing_angle_in)
        y2_oe = numpy.zeros_like(x2_oe)
    else:
        a = numpy.loadtxt(error_file)
        x2_oe = a[:, 0]
        y2_oe = a[:, 1]

    if shape == 0:
        pass
    elif shape == 1:
        if radius >= 0:
            height = radius - numpy.sqrt(radius ** 2 - x2_oe ** 2)
        else:
            height = radius + numpy.sqrt(radius ** 2 - x2_oe ** 2)
        y2_oe += height
    else:
        raise Exception("Wrong shape")

    output_wavefront = propagator1D_offaxis(input_wavefront, x2_oe, y2_oe,
                                            p_distance, q_distance,
                                            grazing_angle_in, grazing_angle_out,
                                            zoom_factor=zoom_factor, normalize_intensities=True)

    # output files
    if write_profile:
        f = open("reflector_profile1D.dat", "w")
        for i in range(height.size):
            f.write("%g %g\n" % (x2_oe[i], y2_oe[i]))
        f.close()
        print("File reflector_profile1D.dat written to disk.")

    return output_wavefront, x2_oe, y2_oe




def calculate_output_wavefront_after_reflector1D(input_wavefront, shape=1, radius=10000.0, grazing_angle=1.5e-3,
                                                 error_flag=0, error_file="", error_edge_management=0, write_profile=0):
    import numpy
    from scipy import interpolate

    output_wavefront = input_wavefront.duplicate()
    abscissas = output_wavefront.get_abscissas()
    abscissas_on_mirror = abscissas / numpy.sin(grazing_angle)

    if shape == 0:
        height = numpy.zeros_like(abscissas_on_mirror)
    elif shape == 1:
        if radius >= 0:
            height = radius - numpy.sqrt(radius ** 2 - abscissas_on_mirror ** 2)
        else:
            height = radius + numpy.sqrt(radius ** 2 - abscissas_on_mirror ** 2)
    else:
        raise Exception("Wrong shape")

    if error_flag:
        a = numpy.loadtxt(error_file)  # extrapolation
        if error_edge_management == 0:
            finterpolate = interpolate.interp1d(a[:, 0], a[:, 1],
                                                fill_value="extrapolate")  # fill_value=(0,0),bounds_error=False)
        elif error_edge_management == 1:
            finterpolate = interpolate.interp1d(a[:, 0], a[:, 1], fill_value=(0, 0), bounds_error=False)
        else:  # crop
            raise Exception("Bad value of error_edge_management")
        height_interpolated = finterpolate(abscissas_on_mirror)
        height += height_interpolated

    phi = -2 * output_wavefront.get_wavenumber() * height * numpy.sin(grazing_angle)

    output_wavefront.add_phase_shifts(phi)

    if error_flag:
        profile_limits = a[-1, 0] - a[0, 0]
        profile_limits_projected = (a[-1, 0] - a[0, 0]) * numpy.sin(grazing_angle)
        wavefront_dimension = output_wavefront.get_abscissas()[-1] - output_wavefront.get_abscissas()[0]
        print("profile deformation dimension: %f m" % (profile_limits))
        print("profile deformation projected perpendicular to optical axis: %f um" % (1e6 * profile_limits_projected))
        print("wavefront window dimension: %f um" % (1e6 * wavefront_dimension))

        if wavefront_dimension <= profile_limits_projected:
            print("\nWavefront window inside error profile domain: no action needed")
        else:
            if error_edge_management == 0:
                print("\nProfile deformation extrapolated to fit wavefront dimensions")
            else:
                output_wavefront.clip(a[0, 0] * numpy.sin(grazing_angle), a[-1, 0] * numpy.sin(grazing_angle))
                print("\nWavefront clipped to projected limits of profile deformation")

    # output files
    if write_profile:
        f = open("reflector_profile1D.dat", "w")
        for i in range(height.size):
            f.write("%g %g\n" % (abscissas_on_mirror[i], height[i]))
        f.close()
        print("File reflector_profile1D.dat written to disk.")

    return output_wavefront, abscissas_on_mirror, height



######################
import time
t0 = time.time()
output_wavefront = calculate_wavefront1D(wavelength=4.9593679373280105e-09,
         wavefront_position=1,
         undulator_length=3.98,
         undulator_distance=13.73,
         x_min=-0.003,
         x_max=0.003,
         number_of_points=3001,
         add_random_phase=False)






input_wavefront = output_wavefront
calculate_output_wavefront_after_reflector1D(input_wavefront,shape=0,
                                             radius=1e+17,
                                             grazing_angle=0.0218,
                                             error_flag=0,
                                             error_file="/home/manuel/Oasys/dabam_profile_140481269429160.dat",
                                             error_edge_management=0,
                                             write_profile=0)

t1 = time.time()

input_wavefront = output_wavefront

# output_wavefront, abscissas_on_mirror, height = calculate_output_wavefront_after_grazing_reflector1D(input_wavefront,
#              radius=1000.0,
#              grazing_angle_in=0.0481860099,
#              grazing_angle_out=0.0727970024,
#              p_distance=12.0,
#              q_distance=4.239,
#              zoom_factor=0.01,
#              error_flag=1,
#              error_file="C:/Users/Manuel/Oasys/VLS_FLEXON.txt",
#              write_profile=0)
output_wavefront, abscissas_on_mirror, height = calculate_output_wavefront_after_grazing_reflector1D(input_wavefront,radius=1000.0,grazing_angle_in=0.0481860099,grazing_angle_out=0.0727970024,p_distance=12.0,q_distance=4.239,zoom_factor=0.01,error_flag=1,error_file="C:/Users/Manuel/Oasys/VLS_FLEXON.txt",write_profile=0)


t2 = time.time()

print("Time 0: ",t1-t0)
print("Time 1: ",t2-t1)
from srxraylib.plot.gol import plot
plot(1e6*output_wavefront.get_abscissas(), output_wavefront.get_intensity())
