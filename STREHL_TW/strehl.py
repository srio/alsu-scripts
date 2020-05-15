import numpy
from srxraylib.plot.gol import plot, set_qt
set_qt()
import scipy.constants as codata

rms = numpy.linspace(0,20,21)
theta = 1.25 * numpy.pi / 180


wavelength = codata.h * codata.c / codata.e / 250
print("Wavelength: %g A" % (wavelength * 1e10))

delta_phi = 2 * rms * 1e-9 * numpy.sin(theta) # wavelength / 14

sr1 = 1 -  (2 * numpy.pi / wavelength * delta_phi)**2
sr2 = numpy.exp(-(2 * numpy.pi / wavelength * delta_phi)**2)
sr3 = (1 -  0.5 * (2 *numpy.pi / wavelength * delta_phi)**2) ** 2
sr4 =  1 - (2 * numpy.pi * rms)**2

print(sr1,sr2,sr3)

plot(rms,sr1,
     rms, sr2,
     rms, sr3,
     legend=["orig","exp","square"],
     yrange=[0,1.1])

