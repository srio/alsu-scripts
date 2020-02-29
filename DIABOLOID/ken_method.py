


import numpy
from oasys.util.oasys_util import write_surface_file


def write_shadow_surface(s,xx,yy,filename='presurface.dat'):
    """
      write_shadowSurface: writes a mesh in the SHADOW/presurface format
      SYNTAX:
           out = write_shadowSurface(z,x,y,filename=filename)
      INPUTS:
           z - 2D array of heights z(x,y)
           x - 1D array of spatial coordinates along mirror width.
           y - 1D array of spatial coordinates along mirror length.

      OUTPUTS:
           filename - output file in SHADOW format. If undefined, the
                     file is names "presurface.dat"

    """

    try:
       fs = open(filename, 'w')
    except IOError:
       out = 0
       print ("Error: can\'t open file: "+filename)
       return
    else:
        # dimensions
        fs.write( "%d  %d \n"%(xx.size,yy.size))
        # y array
        for i in range(yy.size):
            fs.write("%20.18g  "%(yy[i]))
        fs.write("\n")
        # for each x element, the x value followed by the corresponding z(y) profile
        for i in range(xx.size):
            tmps = ""
            for j in range(yy.size):
                tmps += "%20.18g  "%(s[i,j])
            fs.write("%20.18g    %s \n"%(xx[i],tmps))
        fs.close()
        print ("write_shadow_surface: File for SHADOW "+filename+" written to disk.")


def ken_diaboloid_point_to_segment(
        p = 29.3,
        q = 19.53,
        theta = 4.5e-3,
        x = numpy.linspace(-0.01,0.01, 101),
        y= numpy.linspace(-0.1, 0.1, 1001),
        detrend=False,
        filename_shadow="",
        filename_h5=""):

    X = numpy.outer(x,numpy.ones_like(y))
    Y = numpy.outer(numpy.ones_like(x),y)

    s = p * numpy.cos(2*theta)
    z0 = p * numpy.sin(2*theta)
    c = p + q

    print("s: %f, z0: %f, c: %f"%(s,z0,c))

    Z = - numpy.sqrt( c**2 + q**2 - s**2 - 2 * Y * (s + q) - 2 * c * numpy.sqrt( X**2 + (q - Y)**2 ) )

    Z += z0
    print(Z.shape,Z.min(),Z.max())

    if detrend:
        zfit = -theta * y
        for i in range(Z.shape[0]):
            Z[i,:] = Z[i,:] - zfit

    # write_surface_file(Z.T, x, y, "C:\\Users\\manuel\\Oasys\\ken_method.h5", overwrite=True)
    if filename_shadow != "":
        write_shadow_surface(Z, x, y, filename=filename_shadow)
        print("SHADOW file %s written to disk."%filename_shadow)

    if filename_h5 != "":
        write_surface_file(Z.T, x, y, filename_h5, overwrite=True)
        print("HDF5 file %s written to disk."%filename_h5)

    return Z, X, Y

def ken_diaboloid_segment_to_point(
        p = 29.3,
        q = 19.53,
        theta = 4.5e-3,
        x = numpy.linspace(-0.01,0.01, 101),
        y= numpy.linspace(-0.1, 0.1, 1001),
        detrend=False,
        filename_shadow="",
        filename_h5=""):

    X = numpy.outer(x,numpy.ones_like(y))
    Y = numpy.outer(numpy.ones_like(x),y)

    s = q * numpy.cos(2*theta)
    z0 = q * numpy.sin(2*theta)
    c = p + q

    print("s: %f, z0: %f, c: %f"%(s,z0,c))

    Z = z0 - numpy.sqrt( 2 * p**2 + z0**2 + 2 * p * q + 2 * (p + s) * Y - 2 * c * numpy.sqrt( X**2 + (Y + p)**2 ) )


    print(Z.shape,Z.min(),Z.max())

    if detrend:
        # zcentral = Z[Z.shape[0]//2,:]
        # zcoeff = numpy.polyfit(y[(y.size//2-10):(y.size//2+10)],
        #                         zcentral[(y.size//2-10):(y.size//2+10)],1)
        # zfit = zcoeff[1] + y * zcoeff[0]
        # print(zcoeff)
        # from srxraylib.plot.gol import plot
        # plot(y,zcentral,y,zfit)
        zfit = theta * y
        for i in range(Z.shape[0]):
            Z[i,:] = Z[i,:] - zfit
        
    # write_surface_file(Z.T, x, y, "C:\\Users\\manuel\\Oasys\\ken_method.h5", overwrite=True)
    if filename_shadow != "":
        write_shadow_surface(Z, x, y, filename=filename_shadow)
        print("SHADOW file %s written to disk."%filename_shadow)

    if filename_h5 != "":
        write_surface_file(Z.T, x, y, filename_h5, overwrite=True)
        print("HDF5 file %s written to disk."%filename_h5)

    return Z, X, Y



if __name__ == "__main__":

    path = "/home/manuel/Oasys/"
    # path = "c:\\Users\\manuel\\Oasys/"

    system = 'wayne'
    # system = 'point_to_segment'
    system = 'bl12.2.2'
    system = 'bl12.2.2_high_demag'

    if system == 'point_to_segment':
        #
        # point to segment
        #
        p = 29.3
        q = 19.53
        theta = 4.5e-3
        y = numpy.linspace(-0.1, 0.1, 1001)
        x = numpy.linspace(-0.01,0.01, 101)

        Z, X, Y = ken_diaboloid_point_to_segment(p,q,theta,x,y,
                                                 detrend=True,
                                                 filename_shadow='',
                                                 filename_h5=path+"diaboloid_ken_point_to_segment.h5")



        from srxraylib.plot.gol import plot_image, plot
        plot_image(Z,1e3*x,1e3*y,aspect='auto')

    elif system == 'wayne':

        # write_surface_file(Z.T, x, y, path_windows+"tmp_ken.h5", overwrite=True)
        #
        # Wayne system
        #

        #
        # point to segment
        #
        source_diaboloid = 19.54
        diaboloid_image = 9.77
        theta = 4.5e-3
        y = numpy.linspace(-0.361, 0.361, 1001)
        x = numpy.linspace(-0.015,0.015, 101)

        Z, X, Y = ken_diaboloid_segment_to_point(source_diaboloid,diaboloid_image,theta,x,y,
                                                 detrend=True,
                                                 filename_shadow='',
                                                 filename_h5=path+"diaboloid_mckinney_ken_method.h5")

        # Z, X, Y = ken_diaboloid_point_to_segment(diaboloid_image,source_diaboloid,theta,x,y,
        #                                          detrend=True,
        #                                          filename_shadow='',
        #                                          filename_h5=path_windows+"diaboloid_mckinney_ken_method.h5")

    elif system == 'bl12.2.2':

        # write_surface_file(Z.T, x, y, path_windows+"tmp_ken.h5", overwrite=True)
        #
        # Wayne system
        #

        #
        # point to segment
        #
        source_diaboloid = 18.8
        diaboloid_image = 26.875-source_diaboloid
        theta = 2e-3
        y = numpy.linspace(-0.4, 0.4, 1001)
        x = numpy.linspace(-0.015,0.015, 101)

        Z, X, Y = ken_diaboloid_segment_to_point(source_diaboloid,diaboloid_image,theta,x,y,
                                                 detrend=True,
                                                 filename_shadow='',
                                                 filename_h5=path+"diaboloid_bl1222_ken_method.h5")

    elif system == 'bl12.2.2_high_demag':

        # write_surface_file(Z.T, x, y, path_windows+"tmp_ken.h5", overwrite=True)
        #
        # Wayne system
        #

        #
        # point to segment
        #
        source_diaboloid = 18.8
        diaboloid_image = 3.76
        theta = 2e-3
        y = numpy.linspace(-0.4, 0.4, 1001)
        x = numpy.linspace(-0.01, 0.01, 101)

        Z, X, Y = ken_diaboloid_segment_to_point(source_diaboloid, diaboloid_image, theta, x, y,
                                                 detrend=True,
                                                 filename_shadow='',
                                                 filename_h5=path + "diaboloid_bl1222_high_demag_ken_method.h5")