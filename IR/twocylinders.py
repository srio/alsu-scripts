

import numpy
from srxraylib.plot.gol import plot, plot_image

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
            fs.write("%g  "%(yy[i]))
        fs.write("\n")
        # for each x element, the x value followed by the corresponding z(y) profile
        for i in range(xx.size):
            tmps = ""
            for j in range(yy.size):
                tmps += "%g  "%(s[i,j])
            fs.write("%g    %s \n"%(xx[i],tmps))
        fs.close()
        print ("write_shadow_surface: File for SHADOW "+filename+" written to disk.")



def twocylinders(y, delta=0.1,radius1=8.84,radius2=8.48,do_plot=False):

    zc1 = numpy.sqrt(radius1**2 - delta**2)
    zc2 = numpy.sqrt(radius2**2 - delta**2)

    A1 = 1.0
    B1 = -2 * zc1
    C1 = y**2 + zc1**2 - radius1**2
    Dis1 = B1**2 - 4*A1*C1
    z1 = (0.5 / A1) * (-B1 - numpy.sqrt(Dis1))

    A2 = 1.0
    B2 = -2 * zc2
    C2 = y**2 + zc2**2 - radius2**2
    Dis2 = B2**2 - 4*A2*C2
    z2 = (0.5 / A2) * (-B2 - numpy.sqrt(Dis2))

    iflat = numpy.argwhere(numpy.abs(y) < delta)
    z1[iflat] = 0
    z2[iflat] = 0

    z = z1.copy()
    iplus = numpy.argwhere(y >= delta)
    z[iplus] = z2[iplus]

    if do_plot:
        plot(y, z1, y, z2, y, z, legend=["z1", "z2", "z"],title="R1: %f, R2: %f"%(radius1,radius2))
    return z

def shadow3file_twocylinders(x=numpy.linspace(-0.02,0.02,11),y=numpy.linspace(-0.25,0.25,200),delta=0.0100,radius1=8.961770,radius2=8.961770,filename="presurface.dat"):

    z = twocylinders(y, delta=delta,radius1=radius1,radius2=radius2)
    Z = numpy.zeros((x.size,y.size))

    for i in range(x.size):
        Z[i,:] = z

    # plot_image(Z,x,y)

    write_shadow_surface(Z,x,y,filename='presurface.dat')

if __name__ == "__main__":
    shadow3file_twocylinders(numpy.linspace(-0.02,0.02,11),
                             numpy.linspace(-0.25,0.25,200),
                             delta=0.0000, radius1=8.961770, radius2=8.961770, filename="presurface.dat")



