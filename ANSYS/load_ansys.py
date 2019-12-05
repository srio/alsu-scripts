import numpy

from scipy import interpolate
from scipy import spatial

import matplotlib.pylab as plt
import matplotlib as mpl

# from scipy.spatial import Delaunay
# from scipy.interpolate import griddata

import h5py
import numpy

from scipy import interpolate
from scipy import spatial

import matplotlib.pylab as plt
import matplotlib as mpl

# from scipy.spatial import Delaunay
# from scipy.interpolate import griddata

import h5py


def load_ansys_file(filename,nx=300,ny=100,xmin=None,xmax=None,ymin=None,ymax=None,
                    four_quadrants=2,do_plot=False):


    a = numpy.loadtxt(filename,skiprows=0)

    node  = numpy.round(a,10)

    # Coordinates

    Xu = node[:,1]  # X=x in m
    Yu = node[:,2]  # Y=z in m
    Zu = node[:,3]  # Z=uy vertical displacement in m


    X0 = Xu + node[:,4]  # X=x in m
    Y0 = Yu + node[:,5]  # Y=z in m
    Z0 = Zu + node[:,6]  # Z=uy vertical displacement in m

    print("X undeformed limits: ", Xu.min(), Xu.max())
    print("Y undeformed limits: ", Yu.min(), Yu.max())
    print("Z undeformed limits: ", Zu.min(), Zu.max())

    print("X limits: ", X0.min(), X0.max())
    print("Y limits: ", Y0.min(), Y0.max())
    print("Z limits: ", Z0.min(), Z0.max())

    if four_quadrants == 2:
        X1 = numpy.concatenate( (X0,-X0) )
        Y1 = numpy.concatenate( (Y0, Y0) )
        Z1 = numpy.concatenate( (Z0, Z0) )

        X = numpy.concatenate( (X1, X1) )
        Y = numpy.concatenate( (Y1,-Y1) )
        Z = numpy.concatenate( (Z1, Z1) )
    else:
        X = X0
        Y = Y0
        Z = Z0

    #
    # Interpolation
    #
    if xmax is None:
        xmax = X.max() # numpy.abs(X).max()
    if ymax is None:
        ymax = Y.max() # numpy.abs(Y).max()
    if xmin is None:
        xmin = X.min()
    if ymin is None:
        ymin = Y.min()

    # triangulation

    Pi = numpy.array([X, Y]).transpose()
    tri = spatial.Delaunay(Pi)

    if do_plot:
        plt.triplot(X0, Y0 , tri.simplices.copy())
        plt.plot(X0, Y0, "or", label = "Data")
        plt.grid()
        plt.legend()
        plt.xlabel("x")
        plt.ylabel("y")
        plt.show()


    if four_quadrants == 2:
        xlin = numpy.linspace(-xmax,xmax,nx)
        ylin = numpy.linspace(-ymax,ymax,ny)
    else:
        xlin = numpy.linspace(xmin,xmax,nx)
        ylin = numpy.linspace(ymin,ymax,ny)

    XLIN =  numpy.outer(xlin,numpy.ones_like(ylin))
    YLIN =  numpy.outer(numpy.ones_like(xlin),ylin)


    # interpolation
    P = numpy.array([XLIN.flatten(), YLIN.flatten() ]).transpose()

    z = interpolate.griddata(Pi, Z, P, method = "cubic").reshape([nx,ny])

    if do_plot:
        plt.contourf(XLIN, YLIN, z, 50, cmap = mpl.cm.jet)
        plt.colorbar()
        plt.contour(XLIN, YLIN, z, 20, colors = "k")
        #plt.triplot(Xi, Yi , tri.simplices.copy(), color = "k")
        plt.plot(X0, Y0, "or", label = "Data")
        plt.legend()
        plt.grid()
        plt.show()


    if four_quadrants == 1:
        z1 = numpy.vstack((numpy.flip(z,axis=0),z[1:,:]))
        z2 = numpy.hstack((numpy.flip(z1,axis=1),z1[:,1:]))

        xlin2 = numpy.concatenate((-xlin[::-1],xlin[1:]))
        ylin2 = numpy.concatenate((-ylin[::-1],ylin[1:]))

        return z2,xlin2,ylin2
    else:
        return z,xlin,ylin

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


def write_h5_surface(s,xx,yy,filename='presurface.hdf5'):
    subgroup_name = "surface_file"
    file = h5py.File(filename, 'w')
    file[subgroup_name + "/X"] = xx
    file[subgroup_name + "/Y"] = yy
    file[subgroup_name + "/Z"] = s.T
    file.close()
    print("write_h5_surface: File for OASYS " + filename + " written to disk.")



if __name__ == "__main__":
    from srxraylib.plot.gol import plot_image, set_qt

    set_qt()

    filename = "s4.txt"
    z,x,y = load_ansys_file(filename,nx=301,ny=51, do_plot=1, four_quadrants=0)

    plot_image(z,x,y,title=filename+" axes as in ANSYS",
               xtitle="X (%d pixels, max:%f)"%(x.size,x.max()),
               ytitle="Y (%d pixels, max:%f)"%(y.size,y.max()),)

    # write_shadow_surface(z,x,y, filename="s4.dat")
    # Note the axes inversion:  that for shadow y is along the beam
    write_h5_surface(z.T,y,x,filename="/Users/srio/Oasys/s4.h5")


