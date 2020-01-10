import numpy
from srxraylib.plot.gol import plot

#
# create input_wavefront
#
#
def source(photon_energy=250):
    from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D
    input_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-0.00075*2,x_max=0.00075*2,number_of_points=5000)
    input_wavefront.set_photon_energy(photon_energy)
    input_wavefront.set_spherical_wave(radius=13.73,center=0,complex_amplitude=complex(1, 0))
    return input_wavefront

def gaussian_profile(wf_in):
    wf = wf_in.duplicate()
    # print(dir(wf))
    A = wf.get_complex_amplitude()
    und_length = 3.98
    p = 13.73
    wavelength = wf.get_wavelength()
    print("wavelength = ", wavelength)

    sigma_photons = .69 * numpy.sqrt(wavelength / und_length)
    print(A.shape)
    X = wf.get_abscissas()
    print(X.shape)
    print(X[-1])
    sigma = p * sigma_photons
    # sigmax = 2 * X[-1,-1] / 10
    # sigmay = 2 * Y[-1,-1] / 10
    # print("Sigmas: ",sigmax,sigmay)
    print("sigma_photons (intensity): ", sigma)
    print("Diffraction limited size (sigma intensity): ", wavelength / 4 / numpy.pi / sigma)
    print("sigma (intensity): ", sigma)
    sigma *= numpy.sqrt(2)
    print("sigma (amplitude): ", sigma)
    Gx = numpy.exp(-X * X / 2 / sigma ** 2)
    wf.set_complex_amplitude(A * Gx)
    return wf

    print("lambda / (2 pi) / sigma_intensity = ", wavelength / 2 / numpy.pi / sigma)

def apply_error_on_M1(wf_in,focal_x=100.000000):

    from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters
    from syned.beamline.beamline_element import BeamlineElement
    from syned.beamline.element_coordinates import ElementCoordinates
    from wofry.propagator.propagators1D.fresnel_zoom import FresnelZoom1D

    input_wavefront = wf_in.duplicate()

    #
    # info on current oe
    #
    #
    #    -------WOIdealLens---------
    #        focal_x: 100.0 m # Focal length in x [horizontal]
    #        focal_y: None m # Focal length in y [vertical]
    #

    #
    # define current oe
    #
    from wofry.beamline.optical_elements.ideal_elements.lens import WOIdealLens

    optical_element = WOIdealLens(name='', focal_x=focal_x, focal_y=None)

    #
    # propagating (***  ONLY THE ZOOM PROPAGATOR IS IMPLEMENTED ***)
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,
                                       coordinates=ElementCoordinates(p=0.000000, q=0.000000,
                                                                      angle_radial=numpy.radians(0.000000),
                                                                      angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront.duplicate(),
                                                   propagation_elements=propagation_elements)
    # self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 0.050000)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name='FRESNEL_ZOOM_1D')

    return output_wavefront


def propagate_from_M1_to_M3(wf_in):
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


    input_wavefront = wf_in.duplicate()

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

def apply_M3_focusing_and_propagete_to_sample(wf_in,focal_x=2.407000):
    #
    # ===== Example of python code to create propagate current element =====
    #

    #
    # Import section
    #
    from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters
    from syned.beamline.beamline_element import BeamlineElement
    from syned.beamline.element_coordinates import ElementCoordinates
    from wofry.propagator.propagators1D.fresnel_zoom import FresnelZoom1D


    input_wavefront = wf_in.duplicate()

    #
    # info on current oe
    #
    #
    #    -------WOIdealLens---------
    #        focal_x: 2.407 m # Focal length in x [horizontal]
    #        focal_y: None m # Focal length in y [vertical]
    #

    #
    # define current oe
    #
    from wofry.beamline.optical_elements.ideal_elements.lens import WOIdealLens

    optical_element = WOIdealLens(name='', focal_x=focal_x, focal_y=None)

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
    propagation_parameters.set_additional_parameters('magnification_x', 0.050000)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name='FRESNEL_ZOOM_1D')

    return output_wavefront


def plot_wavefront_intensity(wf):
    plot(1e6*wf.get_abscissas(), wf.get_intensity(),
         title="FWHM = %f um, I0=%f"%(1e6*get_wavefront_intensity_fwhm(wf),
                                      get_wavefront_intensity_I0(wf)),
         xtitle="X [um]",ytitle="Intensity [a.u.]")

def get_wavefront_intensity_fwhm(wf):
    from oasys.util.oasys_util import get_fwhm
    fwhm, quote, coordinates = get_fwhm(wf.get_intensity(),wf.get_abscissas())
    return fwhm

def get_wavefront_intensity_I0(wf):
    I = wf.get_intensity()
    return I[I.size//2]


def get_R_incidence_deg(p_foc,q_foc,incidence_deg):
    mm = (1.0 / p_foc + 1.0 / q_foc)
    return 2 / (numpy.cos(incidence_deg * numpy.pi / 180)) / mm

def get_R_grazing(p_foc,q_foc,grazing):
    mm = (1.0 / p_foc + 1.0 / q_foc)
    return 2 / (numpy.sin(grazing)) / mm

def RUN_WOFRY(photon_energy=250,do_plot=True,do_optimize_M3=False,error_radius=1e10):

    wf = gaussian_profile(source(photon_energy=photon_energy))
    # plot(wf.get_abscissas(),wf.get_intensity())


    focal_x_error = 0.5 * error_radius * numpy.sin(1.25 * numpy.pi / 180)
    wf0 = apply_error_on_M1(wf, focal_x=focal_x_error)

    wf1 = propagate_from_M1_to_M3(wf0)
    # plot(wf1.get_abscissas(), wf1.get_intensity())

    #
    # optimization (adaptive optics)
    #
    focal_x = 1.0/(1.0/(13.73+13.599)+1.0/(2.64)) # 2.407000

    if do_optimize_M3:
        DELTA_FOCAL_X = numpy.linspace(-0.4,0.4,200)
        FWHM = numpy.zeros_like(DELTA_FOCAL_X)
        # I0 = numpy.zeros_like(DELTA_FOCAL_X)

        for i in range(DELTA_FOCAL_X.size):
            wf2 = apply_M3_focusing_and_propagete_to_sample(wf1,focal_x=focal_x+DELTA_FOCAL_X[i])
            # plot_wavefront_intensity(wf2)
            FWHM[i] = get_wavefront_intensity_fwhm(wf2)
            # I0[i] = get_wavefront_intensity_I0(wf2)
            # print("FWHM is: ",FWHM[i])

        if do_plot:
            plot(DELTA_FOCAL_X, FWHM, ytitle="FWHM", xtitle="Delta focal distance [m]")
            plot(DELTA_FOCAL_X, I0, ytitle="I0", xtitle="Delta focal distance [m]")

        i_opt = FWHM.argmin()
        print("Start FWHM [um]: %f, optimized: %f : "%(focal_x,focal_x+FWHM[i_opt]))
        # print("Radius [start]: %f",get_R_incidence_deg(13.73+13.599,2.64,90-1.25))
        print("Radius initial: %f optimized: %f "%
              (get_R_incidence_deg(1./(1/focal_x-1./2.64), 2.64, 90 - 1.25),
              get_R_incidence_deg(1./(1/(focal_x+FWHM[i_opt])-1./2.64), 2.64, 90 - 1.25)))

        # calculate optimized
        wf2 = apply_M3_focusing_and_propagete_to_sample(wf1,focal_x=focal_x+DELTA_FOCAL_X[i_opt])
    else:
        wf2 = apply_M3_focusing_and_propagete_to_sample(wf1, focal_x=focal_x)

    if do_plot:
        plot_wavefront_intensity(wf2)

    print("Using error with radius: %f focal: %f"%(error_radius,focal_x_error))

    return wf2

if __name__ == "__main__":
    # wf2 = RUN_WOFRY(photon_energy=250,do_plot=False,do_optimize_M3=False,error_radius=1e6)
    # plot_wavefront_intensity(wf2)

    # wf2 = RUN_WOFRY(photon_energy=250,do_plot=False,do_optimize_M3=True,error_radius=1e2)
    # plot_wavefront_intensity(wf2)


    ERROR_RADIUS = numpy.logspace(1,6,100)
    I0uncorr = numpy.zeros_like(ERROR_RADIUS)
    I0corr = numpy.zeros_like(ERROR_RADIUS)
    I0uncorr1500 = numpy.zeros_like(ERROR_RADIUS)
    I0corr1500 = numpy.zeros_like(ERROR_RADIUS)

    factor = -1.0 # -1.0
    for i in range(ERROR_RADIUS.size):
        wf2 = RUN_WOFRY(photon_energy=250, do_plot=False, do_optimize_M3=True, error_radius=factor*ERROR_RADIUS[i])
        I0corr[i] = get_wavefront_intensity_I0(wf2)
        wf2 = RUN_WOFRY(photon_energy=250, do_plot=False, do_optimize_M3=False, error_radius=factor*ERROR_RADIUS[i])
        I0uncorr[i] = get_wavefront_intensity_I0(wf2)
        wf2 = RUN_WOFRY(photon_energy=1500, do_plot=False, do_optimize_M3=True, error_radius=factor*ERROR_RADIUS[i])
        I0corr1500[i] = get_wavefront_intensity_I0(wf2)
        wf2 = RUN_WOFRY(photon_energy=1500, do_plot=False, do_optimize_M3=False, error_radius=factor*ERROR_RADIUS[i])
        I0uncorr1500[i] = get_wavefront_intensity_I0(wf2)

    plot(ERROR_RADIUS,I0uncorr/I0uncorr[-1],
         ERROR_RADIUS,I0corr/I0corr[-1],
         ERROR_RADIUS, I0uncorr1500 / I0uncorr1500[-1],
         ERROR_RADIUS, I0corr1500 / I0corr1500[-1],
         xlog=True,
         legend=["Uncorrected E=250eV","Corrected E=250 eV","Uncorrected E=1500eV","Corrected E=1500 eV"],
         xtitle="Radius [m]",ytitle="Strehl I/I0")