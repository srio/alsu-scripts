"""

fresnel: 

        functions: 
             goFromTo: calculates the phase shift matrix
 

"""

__author__ = "Manuel Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright = "ESRF, 2012"

import numpy

def goFromToMatrix(field1, source, image, distance=1.0, wavelength=1e-10):
    distance = numpy.array(distance)
    x1 = numpy.outer(source,numpy.ones(image.size))
    x2 = numpy.outer(numpy.ones(source.size),image)
    r = numpy.sqrt( numpy.power(x1-x2,2) + numpy.power(distance,2) )
    wavenumber = numpy.pi*2/wavelength
    T = numpy.exp(1.j * wavenumber *  r)
    return numpy.dot(field1, T)


def goFromToSequential(field1,x1,y1,x2,y2,wavelength=1e-10,normalize_intensities=False):
    field2 = numpy.zeros_like(x2,dtype=complex)
    wavenumber = numpy.pi * 2 / wavelength

    for i in range(field2.size):
        r = numpy.sqrt( numpy.power(x1-x2[i],2) + numpy.power(y1-y2[i],2) )
        field2[i] = ( field1 * numpy.exp(1.j * wavenumber * r) ).sum()

    if normalize_intensities:
        field2 *= numpy.sqrt( (numpy.abs(field1)**2).sum() / (numpy.abs(field2)**2).sum() )
    return field2


def propagator1D_offaxis(input_wavefront,x2_oe,y2_oe,p,q,theta_grazing_in,theta_grazing_out=None,
                         zoom_factor=1.0,normalize_intensities=False):

    from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D

    if theta_grazing_out is None:
        theta_grazing_out = theta_grazing_in

    x1 = input_wavefront.get_abscissas()
    field1 = input_wavefront.get_complex_amplitude()
    wavelength = input_wavefront.get_wavelength()

    x1_oe = -p * numpy.cos(theta_grazing_in) + x1 * numpy.sin(theta_grazing_in)
    y1_oe = p * numpy.sin(theta_grazing_in) + x1 * numpy.cos(theta_grazing_in)

    print(a.shape)

    # field2 is the electric field in the mirror
    field2 = goFromToSequential(field1, x1_oe, y1_oe, x2_oe, y2_oe,
                                wavelength=wavelength,normalize_intensities=normalize_intensities)

    x3 = x1 * zoom_factor

    x3_oe = q * numpy.cos(theta_grazing_out) + x3 * numpy.sin(theta_grazing_out)
    y3_oe = q * numpy.sin(theta_grazing_out) + x3 * numpy.cos(theta_grazing_out)

    # field2 is the electric field in the image plane
    field3 = goFromToSequential(field2, x2_oe, y2_oe, x3_oe, y3_oe,
                                wavelength=wavelength, normalize_intensities=normalize_intensities)


    output_wavefront = GenericWavefront1D.initialize_wavefront_from_arrays(x3, field3, wavelength=wavelength)

    return output_wavefront


def diffraction_by_aperture():

    wavelength = 1.24e-10 # 10keV
    aperture_diameter = 40e-6 # 1e-3 # 1e-6
    detector_size = 300e-6
    distance = 3.6


    sourcepoints = 1000
    detpoints =  2000

    sourcesize = aperture_diameter
    
    position1x = numpy.linspace(-sourcesize/2,sourcesize/2,sourcepoints)
    position2x = numpy.linspace(-detector_size/2,detector_size/2,detpoints)
    field1 = numpy.ones_like(position1x,dtype=complex)
    
    # fieldComplexAmplitude = goFromToMatrix(field1, position1x, position2x, distance, wavelength=wavelength)
    field2 = goFromToSequential(field1,
                                position1x, position1x * 0,
                                position2x, position2x * 0 + distance, wavelength=wavelength)

    #prepare results
    print ("Shape of Complex U2: ", field2.shape)
    print ("Shape of position1x: ",position1x.shape)
    fieldIntensity = numpy.power(numpy.abs(field2), 2)
    fieldPhase = numpy.arctan2(numpy.real(field2), \
                               numpy.imag(field2))

    #
    #plots
    #
    from matplotlib import pylab as plt

    plt.figure(1)
    plt.plot(position2x*1e6,fieldIntensity)
    plt.title("Fresnel-Kirchhoff Diffraction")
    plt.xlabel("X [um]")
    plt.ylabel("Intensity [a.u.]")
    plt.show()


if __name__ == '__main__':
    # diffraction_by_aperture()

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


    output_wavefront = calculate_wavefront1D(wavelength=4.9593679373280105e-09,
                                             wavefront_position=1,
                                             undulator_length=3.98,
                                             undulator_distance=13.73,
                                             x_min=-0.00147,
                                             x_max=0.00147,
                                             number_of_points=3001,
                                             add_random_phase=0)
    from srxraylib.plot.gol import plot

    # plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity())

    #
    # propagator
    #

    input_wavefront = output_wavefront.duplicate()

    p = 1.25
    q = 15.0
    theta_grazing_in = 0.0015
    # theta_grazing_out = 0.0015
    zoom_factor = 0.5

    profile_file = "/home/manuel/Oasys/correction_profile1D.dat"  # "C:\\Users\\manuel\\Oasys\\correction_profile1D.dat"
    a = numpy.loadtxt(profile_file)
    x2_oe = a[:,0]
    y2_oe = a[:,1]

    output_wavefront = propagator1D_offaxis(input_wavefront, x2_oe, y2_oe, p, q, theta_grazing_in, theta_grazing_out=None,
                         zoom_factor=zoom_factor,normalize_intensities=True)


    # plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity())

    print(a.shape)
    print("wf1 intensity", input_wavefront.get_intensity().sum())
    print("wf2 intensity", output_wavefront.get_intensity().sum())
