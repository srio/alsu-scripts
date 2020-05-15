import numpy
from shadow4.optical_surfaces.conic import Conic
from numba import jit, prange

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

def propagate_in_vacuum(input_wavefront,magnification_x=10.0,propagator_flag='z'):
    #
    # Import section
    #
    import numpy
    from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters
    from syned.beamline.beamline_element import BeamlineElement
    from syned.beamline.element_coordinates import ElementCoordinates
    from wofry.propagator.propagators1D.fresnel_zoom import FresnelZoom1D
    from wofry.propagator.propagators1D.integral import Integral1D
    #
    # info on current oe
    #
    #
    #    -------WOScreen1D---------
    #        -------BoundaryShape---------
    #

    #
    # define current oe
    #
    from wofry.beamline.optical_elements.ideal_elements.screen import WOScreen1D

    optical_element = WOScreen1D()

    #
    # propagating (***  ONLY THE ZOOM PROPAGATOR IS IMPLEMENTED ***)
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,
                                       coordinates=ElementCoordinates(p=15.000000, q=0.000000,
                                                                      angle_radial=numpy.radians(0.000000),
                                                                      angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront.duplicate(),
                                                   propagation_elements=propagation_elements)
    # self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', magnification_x)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
        propagator.add_propagator(Integral1D())
    except:
        pass

    if propagator_flag == 'z':
        handler_name = 'FRESNEL_ZOOM_1D'
    elif propagator_flag == 'i':
        handler_name = 'INTEGRAL_1D'

    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name=handler_name)

    return output_wavefront


def calculate_output_wavefront_after_reflector1D(input_wavefront, shape=1, radius=10000.0, wavy_amplitude=1e-9,
                                                 wavy_ripples=1.0, grazing_angle=1.5e-3, error_flag=0, error_file="",
                                                 error_edge_management=0, write_profile=""):
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
    elif shape > 1:
        if wavy_ripples == 0.0:
            height = numpy.zeros_like(abscissas_on_mirror)
        else:
            period = (abscissas_on_mirror[-1] - abscissas_on_mirror[0]) / wavy_ripples
            if shape == 2:
                y = numpy.cos(2 * numpy.pi * abscissas_on_mirror / period)
            elif shape == 3:
                y = numpy.sin(2 * numpy.pi * abscissas_on_mirror / period)
            y -= y.min()
            y /= y.max()
            y *= wavy_amplitude
            height = y
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
    if write_profile != "":
        f = open(write_profile, "w")
        for i in range(height.size):
            f.write("%g %g\n" % (abscissas_on_mirror[i], height[i]))
        f.close()
        print("File %s written to disk." % write_profile)

    return output_wavefront, abscissas_on_mirror, height


def run_beamline(error_flag=1,error_file="/Users/srio/Oasys/dabam_profile_6010353992.dat"):


    # wf0 = calculate_wavefront1D(wavelength=4.9593679373280105e-09,
    #                                          wavefront_position=0,
    #                                          undulator_length=2.24,
    #                                          undulator_distance=13.73,
    #                                          x_min=-0.00025,
    #                                          x_max=0.00025,
    #                                          number_of_points=1000,
    #                                          add_random_phase=0)

    #
    # create input_wavefront
    #
    #
    from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D
    input_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-0.0002, x_max=0.0002,
                                                                         number_of_points=1000)
    input_wavefront.set_photon_energy(250)
    input_wavefront.set_gaussian(sigma_x=0.000012, amplitude=1.000000, shift=0.000000)
    wf0 = input_wavefront

    # from srxraylib.plot.gol import plot
    # plot(input_wavefront.get_abscissas(), input_wavefront.get_intensity())







    #
    #
    #
    wf1 = propagate_in_vacuum(wf0)

    # plot(1e6 * wf1.get_abscissas(), wf1.get_intensity())

    #
    #
    #
    wf2, abscissas_on_mirror, height = calculate_output_wavefront_after_reflector1D(wf1,
            shape=1,radius=687.6038987249137,grazing_angle=0.02181661,
            error_flag=error_flag,
            error_file=error_file,
            error_edge_management=0,write_profile="")


    # plot(1e6 * wf2.get_abscissas(),wf2.get_intensity())


    #
    #
    #
    wf3 = propagate_in_vacuum(wf2,magnification_x=0.03,propagator_flag='i')

    return wf3




@jit(nopython=True, parallel=True)
def goFromToSequential(field1, x1, y1, x2, y2, wavelength=1e-10, normalize_intensities=False):
    field2 = x2 * 0j
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

    output_wavefront = GenericWavefront1D.initialize_wavefront_from_arrays(x3, field3  / numpy.sqrt(zoom_factor), wavelength=wavelength)

    return output_wavefront

def calculate_output_wavefront_after_grazing_reflector1D(input_wavefront,shape=1,
                                        p_focus=1.0,
                                        q_focus=1.0,
                                        grazing_angle_in=1.5e-3,
                                        p_distance=1.0,
                                        q_distance=1.0,
                                        zoom_factor=1.0,
                                        error_flag=0, error_file="", write_profile=0):



    x1 = input_wavefront.get_abscissas()
    field1 = input_wavefront.get_complex_amplitude()

    if error_flag == 0: # no profile file
        x2_oe = x1 / numpy.sin(grazing_angle_in)
        y2_oe = numpy.zeros_like(x2_oe)
    else:
        a = numpy.loadtxt(error_file)
        x2_oe = a[:, 0]
        y2_oe = a[:, 1]


    if shape == 0:
        pass
    elif shape == 1:
        ccc = Conic.initialize_as_sphere_from_focal_distances(p_focus, q_focus, grazing_angle_in)
        height = ccc.height(x2_oe)
        print(ccc.info())
        y2_oe += height
    elif shape == 2:
        ccc = Conic.initialize_as_ellipsoid_from_focal_distances(p_focus, q_focus, grazing_angle_in)
        height = ccc.height(x2_oe)
        print(ccc.info())
        y2_oe += height
    elif shape == 3:
        ccc = Conic.initialize_as_paraboloid_from_focal_distances(p_focus, q_focus, grazing_angle_in)
        height = ccc.height(x2_oe)
        print(ccc.info())
        y2_oe += height
    else:
        raise Exception("Wrong shape")

    output_wavefront = propagator1D_offaxis(input_wavefront, x2_oe, y2_oe,
                                                 p_distance,q_distance,
                                                 grazing_angle_in,
                                                 zoom_factor=zoom_factor,normalize_intensities=True)

    # output files
    if write_profile:
        f = open("reflector_profile1D.dat","w")
        for i in range(height.size):
            f.write("%g %g\n"%(x2_oe[i],y2_oe[i]))
        f.close()
        print("File reflector_profile1D.dat written to disk.")

    return output_wavefront, x2_oe, y2_oe

def run_beamline_2(error_flag=0,error_file=""):

    #
    # main

    # create input_wavefront
    #
    #
    from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D
    input_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-0.0002, x_max=0.0002,
                                                                         number_of_points=5000)
    input_wavefront.set_photon_energy(250)
    input_wavefront.set_gaussian(sigma_x=0.000012, amplitude=1.000000, shift=0.000000)


    output_wavefront, abscissas_on_mirror, height = calculate_output_wavefront_after_grazing_reflector1D(
        input_wavefront, shape=2, p_focus=15.0, q_focus=15.0, grazing_angle_in=0.02181661, p_distance=15.0,
        q_distance=15.0,
        zoom_factor=10, #0.3,
        error_flag=error_flag,
        error_file=error_file, write_profile=0)

    # from srxraylib.plot.gol import plot
    # plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity())
    return output_wavefront

if __name__ == "__main__":
    from srxraylib.plot.gol import plot, plot_image, set_qt
    set_qt()


    wf = run_beamline(error_flag=1,error_file="C:\\Users\\Manuel\\Oasys/dabam_profile_1728934226968.dat")
    wf2 = run_beamline_2(error_flag=1, error_file="C:\\Users\\Manuel\\Oasys/dabam_profile_1728934226968.dat")

    import matplotlib
    plot(-1e3 * wf.get_abscissas(), wf.get_intensity(),
         1e3 * wf2.get_abscissas(), wf2.get_intensity(),
         show=0)
    matplotlib.pyplot.grid()
    matplotlib.pyplot.show()