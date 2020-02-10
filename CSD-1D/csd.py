
import numpy
from srxraylib.plot.gol import plot_image, plot
from wofry.propagator.wavefront2D.generic_wavefront import GenericWavefront2D


def get_fwhm(histogram, bins):
    quote = numpy.max(histogram)*0.5
    cursor = numpy.where(histogram >= quote)

    if histogram[cursor].size > 1:
        bin_size    = bins[1]-bins[0]
        fwhm        = bin_size*(cursor[0][-1]-cursor[0][0])
        coordinates = (bins[cursor[0][0]], bins[cursor[0][-1]])
    else:
        fwhm = 0.0
        coordinates = None

    return fwhm, quote, coordinates


def csd(X1,X2,sigmaI,sigmaMu):

    return numpy.exp( - (X1**2 + X2**2) / (4*sigmaI**2)) * numpy.exp(- (X1 - X2)**2 / 2 /sigmaMu**2) + 0j


def fresnel(wavefront,propagation_distance,shift_half_pixel=False):
    wavelength = wavefront.get_wavelength()

    #
    # convolving with the Fresnel kernel via FFT multiplication
    #
    fft = numpy.fft.fft2(wavefront.get_complex_amplitude())

    # frequency for axis 1
    shape = wavefront.size()
    delta = wavefront.delta()

    pixelsize = delta[0] # p_x[1] - p_x[0]
    npixels = shape[0]
    freq_nyquist = 0.5/pixelsize
    freq_n = numpy.linspace(-1.0,1.0,npixels)
    freq_x = freq_n * freq_nyquist

    # frequency for axis 2
    pixelsize = delta[1]
    npixels = shape[1]
    freq_nyquist = 0.5/pixelsize
    freq_n = numpy.linspace(-1.0,1.0,npixels)
    freq_y = freq_n * freq_nyquist

    if shift_half_pixel:
        freq_x = freq_x - 0.5 * numpy.abs(freq_x[1] - freq_x[0])
        freq_y = freq_y - 0.5 * numpy.abs(freq_y[1] - freq_y[0])

    freq_xy = numpy.array(numpy.meshgrid(freq_y,freq_x))
    # fft *= numpy.exp((-1.0j) * numpy.pi * wavelength * propagation_distance *
    #               numpy.fft.fftshift(freq_xy[0]*freq_xy[0] + freq_xy[1]*freq_xy[1]) )
    fft *= numpy.exp((-1.0j) * numpy.pi * wavelength * propagation_distance *
                  numpy.fft.fftshift(freq_xy[0]*freq_xy[0] + freq_xy[1]*freq_xy[1]) )

    wf_propagated = GenericWavefront2D.initialize_wavefront_from_arrays(x_array=wavefront.get_coordinate_x(),
                                                                        y_array=wavefront.get_coordinate_y(),
                                                                        z_array=numpy.fft.ifft2(fft),
                                                                        wavelength=wavelength)

    return wf_propagated

def fresnel_zoom(wavefront1,propagation_distance,magnification_x=1.0,magnification_y=1.0,shift_half_pixel=False):


    wavefront = wavefront1.duplicate()
    wavelength = wavefront.get_wavelength()
    wavenumber = wavefront.get_wavenumber()

    shape = wavefront.size()
    delta = wavefront.delta()

    pixelsize = delta[0]
    npixels = shape[0]
    freq_nyquist = 0.5 / pixelsize
    freq_n = numpy.linspace(-1.0, 1.0, npixels)
    freq_x = freq_n * freq_nyquist

    # frequency for axis 2
    pixelsize = delta[1]
    npixels = shape[1]
    freq_nyquist = 0.5 / pixelsize
    freq_n = numpy.linspace(-1.0, 1.0, npixels)
    freq_y = freq_n * freq_nyquist

    if shift_half_pixel:
        freq_x = freq_x - 0.5 * numpy.abs(freq_x[1] - freq_x[0])
        freq_y = freq_y - 0.5 * numpy.abs(freq_y[1] - freq_y[0])

    f_x, f_y = numpy.meshgrid(freq_x, freq_y, indexing='ij')
    # fsq = numpy.fft.fftshift(f_x ** 2 / magnification_x + f_y ** 2 / magnification_y)
    fsq = numpy.fft.fftshift(f_x ** 2 / magnification_x + f_y ** 2 / magnification_y)

    x = wavefront.get_mesh_x()
    y = wavefront.get_mesh_y()

    x_rescaling = wavefront.get_mesh_x() * magnification_x
    y_rescaling = wavefront.get_mesh_y() * magnification_y

    r1sq = x ** 2 * (1 - magnification_x) + y ** 2 * (1 - magnification_y)
    r2sq = x_rescaling ** 2 * ((magnification_x - 1) / magnification_x) + y_rescaling ** 2 * ((magnification_y - 1) / magnification_y)

    Q1 = wavenumber / 2 / propagation_distance * r1sq
    Q2 = numpy.exp(-1.0j * numpy.pi * wavelength * propagation_distance * fsq)
    Q3 = numpy.exp(1.0j * wavenumber / 2 / propagation_distance * r2sq)

    wavefront.add_phase_shift(Q1)

    fft = numpy.fft.fft2(wavefront.get_complex_amplitude())

    ifft = numpy.fft.ifft2(fft * Q2) * Q3 / numpy.sqrt(magnification_x * magnification_y)

    wf_propagated = GenericWavefront2D.initialize_wavefront_from_arrays(x_array=wavefront.get_coordinate_x()*magnification_x,
                                                                        y_array=wavefront.get_coordinate_y()*magnification_y,
                                                                        z_array=ifft,
                                                                        wavelength=wavelength)
    return wf_propagated

def get_intensity_from_csd(csd):
    intensity = numpy.zeros(csd.shape[0])
    for i in range(intensity.size):
        intensity[i] = csd[i,i]
    return intensity

if __name__ == "__main__":

    beta = 5
    sigmaI = 1.5e-6
    sigmaMu = beta * sigmaI
    propagation_distance = 10.0
    do_plot = 0


    sigmaLim = 2 * max([sigmaI,sigmaMu])

    x = numpy.linspace(-10e-6, 10e-6, 1500) #-sigmaLim,sigmaLim,1500)

    X1 = numpy.outer(x,numpy.ones_like(x))
    X2 = numpy.outer(numpy.ones_like(x),x)

    # plot_image(X2,x,x)
    z = csd(X1,X2,sigmaI,sigmaMu)
    zi = get_intensity_from_csd(z)
    if do_plot:
        plot_image(z, 1e6 * x, 1e6 * x, xtitle="X1", ytitle="X2", title="CDS before propagation", show=0)
        plot(1e6 * x, zi,show=0)

    wf = GenericWavefront2D.initialize_wavefront_from_arrays(x,x,z)
    wf.set_photon_energy(23000)

    # wfp = fresnel(wf,propagation_distance=propagation_distance)
    magnification = 5
    wfp = fresnel_zoom(wf,propagation_distance=propagation_distance,magnification_x=magnification,magnification_y=magnification)


    zp = wfp.get_intensity()
    zpi =  get_intensity_from_csd(zp)

    fwhm_z, tmp, tmp = get_fwhm(zi / zi.max(), 1e6 * x)
    fwhm_zp, tmp, tmp = get_fwhm(zpi / zpi.max(), 1e6 * x)

    if do_plot:
        plot_image(zp, 1e6 * wfp.get_coordinate_x(), 1e6 * wfp.get_coordinate_y(), xtitle="X1", ytitle="X2",
                   title="CSD after propagation", show=0)
        plot(1e6 * x, zi/zi.max(), 1e6 * x, zpi/zpi.max(),
             legend=["before propagation fwhm=%f"%fwhm_z,"after propagation fwhm=%f"%fwhm_zp])

    sigma_theory = wf.get_wavelength() / 4 / numpy.pi / sigmaI * propagation_distance
    print("Theoretical propagated sigma: %f um, calculated: %f um"%(1e6*sigma_theory,fwhm_zp/2.35* numpy.sqrt(1)) )
    print("Theoretical propagated FWHM: %f um, calculated: %f um" % (1e6 * 2.35 * sigma_theory, fwhm_zp * numpy.sqrt(1)))