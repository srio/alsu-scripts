import scipy.constants as codata
import numpy
from srxraylib.plot.gol import plot

energy = 0.01 # 15000.0
wavelength = codata.h * codata.c / energy / codata.e

print("energy: %f eV, wavelength: %f um"%(energy,wavelength*1e6))

k = 2 * numpy.pi / wavelength

x = numpy.linspace(-20e-3,20e-3,500)

divergence = 25e-3 # 1e-3
arg = k * x / 2 * divergence

y = (numpy.sin(arg) / arg )**2


t = numpy.argwhere(y/y.max() > 0.5)
# print(t)
print(t[-1] , t[0])
print(x[[-1]] , x[[0]])

fwhm = x[t[-1]] - x[t[0]]
print("fwhm: %g"%fwhm)



print("fwhm = ", 2 / divergence / k * 2.78)

print("divergence: %g rad, fwhm = %g "%(divergence, 0.8849 * wavelength / divergence))
print("Rayleigh: %g rad, fwhm = %g "%(divergence, 0.61 * wavelength / divergence))

divergence = 0.035
print("\ndivergence: %g rad, fwhm = %g "%(divergence, 0.8849 * wavelength / divergence))
print("Rayleigh: %g rad, fwhm = %g "%(divergence, 0.61 * wavelength / divergence))

