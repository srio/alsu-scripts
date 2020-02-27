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


def goFromToSequential(field1,x1,y1,x2,y2,wavelength=1e-10):
    field2 = numpy.zeros_like(x2,dtype=complex)
    wavenumber = numpy.pi * 2 / wavelength

    for i in range(field2.size):
        r = numpy.sqrt( numpy.power(x1-x2[i],2) + numpy.power(y1-y2[i],2) )
        field2[i] = ( field1 * numpy.exp(1.j * wavenumber * r) ).sum()

    return field2


if __name__ == '__main__':

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
    # write spec formatted file
    #
    # out_file = "fresnel_kirchhoff_1D.spec"
    # f = open(out_file, 'w')
    # header="#F %s \n\n#S  1 fresnel-kirchhoff diffraction integral\n#N 3 \n#L X[m]  intensity  phase\n"%out_file
    #
    # f.write(header)
    #
    # for i in range(detpoints):
    #    out = numpy.array((position2x[i], fieldIntensity[i], fieldPhase[i]))
    #    f.write( ("%20.11e "*out.size+"\n") % tuple( out.tolist())  )
    #
    # f.close()
    # print ("File written to disk: %s"%out_file)

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