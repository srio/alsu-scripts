




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



def get_wavefront_intensity_fwhm(wf):
    from oasys.util.oasys_util import get_fwhm
    fwhm, quote, coordinates = get_fwhm(wf.get_intensity(),wf.get_abscissas())
    return fwhm

def get_wavefront_intensity_I0(wf):
    I = wf.get_intensity()
    return I[I.size//2]

def do_propagation(input_wavefront,q=2.640000):

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
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=0.000000,    q=q,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront.duplicate(),    propagation_elements = propagation_elements)
    #self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 0.010000)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='FRESNEL_ZOOM_1D')

    return output_wavefront


#
# create/import your input_wavefront
#
#
from srxraylib.plot.gol import plot
from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D
input_wavefront = GenericWavefront1D.initialize_wavefront_from_range(-10e-6,10e-6,200,1e-10)

input_wavefront = GenericWavefront1D.load_h5_file("/home/manuel/Oasys/tmp.h5")
# plot(input_wavefront.get_abscissas(),input_wavefront.get_intensity(),title="input")


q=2.640000
DELTAQ = numpy.linspace(-0.1,0.1,100)
Q = DELTAQ + q
I0 = numpy.zeros_like(DELTAQ)
for i,deltaq in enumerate(DELTAQ):
    output_wavefront = do_propagation(input_wavefront,q=q+deltaq)
    I0[i] = get_wavefront_intensity_I0(output_wavefront)
    # plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title=I0[i])


plot(Q,I0)
imax = I0.argmax()
print(Q[imax],I0[imax])