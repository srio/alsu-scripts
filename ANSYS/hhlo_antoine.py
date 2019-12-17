def ansys_read(filename, nx=101, ny=51):
    """
    Read an High heatload deformation files (ALS-U standards)

    :param filename: full filename of the file to read
    :param nx: number
    :returns: (X, Y, Z) 2D arrays of (X,Y) locations and corresponding height Z, in meter
    """

    import numpy
    import csv
    import scipy.interpolate as interpolate

    # read data wih csv
    with open(filename, 'r') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=" ", skipinitialspace=True)
        row_count = sum(1 for row in csvreader)

    data = numpy.zeros((row_count, 7))

    with open(filename, 'r') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=" ", skipinitialspace=True)
        i = 0
        for row in csvreader:
            data[i, :] = numpy.asarray(row).astype(numpy.float)
            i = i + 1

    # parse data
    x = data[:, 1]
    dx = data[:, 4]
    y = data[:, 2]
    dy = data[:, 5]
    z = data[:, 3]
    dz = data[:, 6]

    # deformation coordinate
    xdx = x + dx
    ydy = y + dy
    zdz = z + dz

    # interpolation on regular grid
    xp = numpy.linspace(numpy.min(x), numpy.max(x), nx)
    yp = numpy.linspace(numpy.min(y), numpy.max(y), ny)
    X, Y = numpy.meshgrid(xp, yp)
    Z = interpolate.griddata((xdx, ydy), zdz, (X, Y), method='linear')

    return X, Y, Z


def ansys_height(filename, nx=101, ny=51):
    import numpy
    Xq, Yq, Zq = ansys_read(filename, nx, ny)
    Nx = Xq.shape[0]
    Ny = Yq.shape[1]
    X_m, Y_m = numpy.meshgrid(numpy.linspace(-Xq[0, -1], Xq[0, -1], 2 * Ny), numpy.linspace(-Yq[-1, 0], Yq[-1, 0], 2 * Nx))
    Z_m = numpy.concatenate((numpy.flip(numpy.concatenate((numpy.flip(Zq, axis=0), Zq), axis=0), axis=1),
                             numpy.concatenate((numpy.flip(Zq, axis=0), Zq), axis=0)), axis=1)
    return X_m, Y_m, Z_m

def write_shadow_surface(s, xx, yy, outFile='presurface.dat'):
    """
      write_shadowSurface: writes a mesh in the SHADOW/presurface format
      SYNTAX:
           out = write_shadowSurface(z,x,y,outFile=outFile)
      INPUTS:
           z - 2D array of heights
           x - 1D array of spatial coordinates along mirror width.
           y - 1D array of spatial coordinates along mirror length.

      OUTPUTS:
           out - 1=Success, 0=Failure
           outFile - output file in SHADOW format. If undefined, the
                     file is names "presurface.dat"

    """
    out = 1

    try:
        fs = open(outFile, 'w')
    except IOError:
        out = 0
        print("Error: can\'t open file: " + outFile)
        return
    else:
        # dimensions
        fs.write(repr(xx.size) + " " + repr(yy.size) + " \n")
        # y array
        for i in range(yy.size):
            fs.write(' ' + repr(yy[i]))
        fs.write("\n")
        # for each x element, the x value and the corresponding z(y)
        # profile
        for i in range(xx.size):
            tmps = ""
            for j in range(yy.size):
                tmps = tmps + "  " + repr(s[j, i])
            fs.write(' ' + repr(xx[i]) + " " + tmps)
            fs.write("\n")
        fs.close()
        print("write_shadow_surface: File for SHADOW " + outFile + " written to disk.")

def compute(filename):
    import numpy as np
    X, Y, Z = ansys_height(filename)
    xx_m = X[1, :]
    yy_m = Y[:, 1]
    write_shadow_surface(np.transpose(Z), yy_m, xx_m, outFile='presurface.dat')
    return np.transpose(Z), yy_m, xx_m


if __name__ == "__main__":
    from srxraylib.plot.gol import plot_image, set_qt
    import os

    set_qt()

    filename = "s4.txt"
    z, x, y = compute(filename)

    plot_image(z,x,y,title=filename+" axes as in ANSYS",
               xtitle="X (%d pixels, max:%f)"%(x.size,x.max()),
               ytitle="Y (%d pixels, max:%f)"%(y.size,y.max()),)
    #
    # # write_shadow_surface(z,x,y, filename="s4.dat")
    # # Note the axes inversion:  that for shadow y is along the beam
    #
    # write_h5_surface(z.T,y,x,filename="/home/manuel/Oasys/s4.h5")