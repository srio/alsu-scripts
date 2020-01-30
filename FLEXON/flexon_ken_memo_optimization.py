import numpy


def get_R_incidence_deg(p_foc,q_foc,incidence_deg):
    mm = (1.0 / p_foc + 1.0 / q_foc)
    return 2 / (numpy.cos(incidence_deg * numpy.pi / 180)) / mm



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


def calculate_output_wavefront_after_reflector1D(input_wavefront, radius=10000.0, grazing_angle=1.5e-3, error_flag=0,
                                                 error_file="", write_profile=0):
    import numpy
    from scipy import interpolate
    output_wavefront = input_wavefront.duplicate()
    abscissas = output_wavefront.get_abscissas()
    abscissas_on_mirror = abscissas / numpy.sin(grazing_angle)
    if radius >= 0:
        height = radius - numpy.sqrt(radius ** 2 - abscissas_on_mirror ** 2)
    else:
        height = radius + numpy.sqrt(radius ** 2 - abscissas_on_mirror ** 2)

    if error_flag:
        a = numpy.loadtxt(error_file)
        finterpolate = interpolate.interp1d(a[:, 0], a[:, 1], fill_value="extrapolate")
        height_interpolated = finterpolate(abscissas_on_mirror)
        height += height_interpolated

    phi = -2 * output_wavefront.get_wavenumber() * height * numpy.sin(grazing_angle)
    output_wavefront.add_phase_shifts(phi)

    # output files
    if write_profile:
        f = open("reflector_profile1D.dat", "w")
        for i in range(height.size):
            f.write("%g %g\n" % (abscissas_on_mirror[i], height[i]))
        f.close()
        print("File reflector_profile1D.dat written to disk.")

    return output_wavefront, abscissas_on_mirror, height


def propagate_to_M3(input_wavefront):
    #
    # ===== Example of python code to create propagate current element =====
    #

    #
    # Import section
    #
    import numpy
    from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters
    from syned.beamline.beamline_element import BeamlineElement
    from syned.beamline.element_coordinates import ElementCoordinates
    from wofry.propagator.propagators1D.fresnel_zoom import FresnelZoom1D


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
                                       coordinates=ElementCoordinates(p=0.000000, q=13.599000,
                                                                      angle_radial=numpy.radians(0.000000),
                                                                      angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront.duplicate(),
                                                   propagation_elements=propagation_elements)
    # self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 2.000000)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name='FRESNEL_ZOOM_1D')

    return output_wavefront


def calculate_output_wavefront_after_reflector1D_M3(input_wavefront, radius=10000.0, grazing_angle=1.5e-3,
                                                 error_flag=0, error_file="", write_profile=0):
    import numpy
    from scipy import interpolate
    output_wavefront = input_wavefront.duplicate()
    abscissas = output_wavefront.get_abscissas()
    abscissas_on_mirror = abscissas / numpy.sin(grazing_angle)
    if radius >= 0:
        height = radius - numpy.sqrt(radius ** 2 - abscissas_on_mirror ** 2)
    else:
        height = radius + numpy.sqrt(radius ** 2 - abscissas_on_mirror ** 2)

    if error_flag:
        a = numpy.loadtxt(error_file)
        finterpolate = interpolate.interp1d(a[:, 0], a[:, 1], fill_value="extrapolate")
        height_interpolated = finterpolate(abscissas_on_mirror)
        height += height_interpolated

    phi = -2 * output_wavefront.get_wavenumber() * height * numpy.sin(grazing_angle)
    output_wavefront.add_phase_shifts(phi)

    # output files
    if write_profile:
        f = open("reflector_profile1D.dat", "w")
        for i in range(height.size):
            f.write("%g %g\n" % (abscissas_on_mirror[i], height[i]))
        f.close()
        print("File reflector_profile1D.dat written to disk.")

    return output_wavefront, abscissas_on_mirror, height


def propagate_to_image(input_wavefront):
    #
    # ===== Example of python code to create propagate current element =====
    #

    #
    # Import section
    #
    import numpy
    from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters
    from syned.beamline.beamline_element import BeamlineElement
    from syned.beamline.element_coordinates import ElementCoordinates
    from wofry.propagator.propagators1D.fresnel_zoom import FresnelZoom1D


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
                                       coordinates=ElementCoordinates(p=0.000000, q=2.640000,
                                                                      angle_radial=numpy.radians(0.000000),
                                                                      angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront.duplicate(),
                                                   propagation_elements=propagation_elements)
    # self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 0.020000)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name='FRESNEL_ZOOM_1D')

    return output_wavefront




if __name__ == "__main__":

    from srxraylib.plot.gol import plot


    #
    # source
    #

    output_wavefront = calculate_wavefront1D(wavelength=4.9593679373280105e-09,
                                             wavefront_position=1,
                                             undulator_length=3.98,
                                             undulator_distance=13.73,
                                             x_min=-0.00147,
                                             x_max=0.00147,
                                             number_of_points=3000,
                                             add_random_phase=False)




    #
    # M1
    #
    input_wavefront = output_wavefront

    output_wavefront, abscissas_on_mirror, height = calculate_output_wavefront_after_reflector1D(input_wavefront,
                                                                                                 radius=496600.0,
                                                                                                 grazing_angle=0.0218,
                                                                                                 error_flag=0,
                                                                                                 error_file="/home/manuel/Oasys/dabam_profile_139713533227224.dat",
                                                                                                 write_profile=1)


    #
    # screen before M3
    #
    input_wavefront = output_wavefront

    output_wavefront = propagate_to_M3(input_wavefront)

#
    #
    # M3 (adaptive)
    #

    x = numpy.linspace(-0.135,0.135,100)
    R = get_R_incidence_deg(13.73+13.599, 2.64, 90 - 1.25)
    print(">>>>>>>>>>>>>>>> R",R,R.shape)
    height = R - numpy.sqrt( R ** 2 - x ** 2)

    f = open("tmp.dat",'w')
    for i in range(x.size):
        f.write("%g %g \n"%(x[i],height[i]))

    f.close()


    input_wavefront = output_wavefront

    output_wavefront, abscissas_on_mirror, height = calculate_output_wavefront_after_reflector1D_M3(input_wavefront,
                                                                                                 radius=220000000000000.72,
                                                                                                 grazing_angle=0.02181,
                                                                                                 error_flag=1,
                                                                                                 error_file="tmp.dat", #"/home/manuel/Oasys/correction_profile1D.dat",
                                                                                                 write_profile=0)

    #
    #
    #

    input_wavefront = output_wavefront

    output_wavefront= propagate_to_image(input_wavefront)


    from srxraylib.plot.gol import plot
    plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity())


