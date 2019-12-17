import numpy

from scipy import interpolate
from scipy import spatial

import matplotlib.pylab as plt
import matplotlib as mpl

def write_generic_h5_surface(s, xx, yy, filename='presurface.hdf5',subgroup_name="surface_file"):
    import h5py
    file = h5py.File(filename, 'w')
    file[subgroup_name + "/X"] = xx
    file[subgroup_name + "/Y"] = yy
    file[subgroup_name + "/Z"] = s.T
    file.close()
    print("write_h5_surface: File for OASYS " + filename + " written to disk.")


class FEA_File():
    def __init__(self,filename=""):
        self.filename = filename
        self.reset()

    def reset(self):
        self.Xundeformed = None  # 1D array
        self.Yundeformed = None  # 1D array
        self.Zundeformed = None  # 1D array

        self.Xdeformation = None  # 1D array
        self.Ydeformation = None  # 1D array
        self.Zdeformation = None  # 1D array

        self.triPi = None
        self.tri = None

        self.x_interpolated = None  # 1D array
        self.y_interpolated = None  # 1D array
        self.Z_INTERPOLATED = None  # 2D array

    @classmethod
    def process_file(cls, filename_in, n_axis_0=301, n_axis_1=51,
                     filename_out="", invert_axes_names=False,
                     detrend=True, reset_height_method=0, do_plot=False):

        o1 = FEA_File(filename=filename_in)
        o1.load_multicolumn_file()

        o1.triangulate()

        if do_plot:
            o1.plot_triangulation()

        o1.interpolate(n_axis_0, n_axis_1)
        if do_plot:
            o1.plot_interpolated()

        if o1.does_interpolated_have_nan():
            o1.remove_borders_in_interpolated_data()

        if do_plot:
            o1.plot_surface_image()

        if detrend:
            o1.detrend()

        # o1.reset_height_to_minimum()

        if reset_height_method == 0:
            pass
        elif reset_height_method == 1:
            o1.reset_height_to_minimum()
        elif reset_height_method == 2:
            o1.reset_height_to_central_value()

        if do_plot:
            o1.plot_surface_image()

        if filename_out != "":
            o1.write_h5_surface(filename=filename_out, invert_axes_names=invert_axes_names)

        return o1

    def set_filename(self,filename):
        self.filename = filename

    def load_multicolumn_file(self,skiprows=0):
        a = numpy.loadtxt(self.filename,skiprows=skiprows )

        node = numpy.round(a, 10)

        # Coordinates

        self.Xundeformed = node[:, 1]  # X=x in m
        self.Yundeformed = node[:, 2]  # Y=z in m
        self.Zundeformed = node[:, 3]  # Z=uy vertical displacement in m

        self.Xdeformation = node[:, 4]  # X=x in m
        self.Ydeformation = node[:, 5]  # Y=z in m
        self.Zdeformation = node[:, 6]  # Z=uy vertical displacement in m

    def Xdeformed(self):
        return self.Xundeformed + self.Xdeformation

    def Ydeformed(self):
        return self.Yundeformed + self.Ydeformation

    def Zdeformed(self):
        return self.Zundeformed + self.Zdeformation


    def get_deformed(self):
        return self.Xdeformed(),self.Ydeformed(),self.Zdeformed()

    def get_limits_undeformed(self):

        print("X undeformed limits: ", self.Xundeformed.min(), self.Xundeformed.max())
        print("Y undeformed limits: ", self.Yundeformed.min(), self.Yundeformed.max())
        print("Z undeformed limits: ", self.Zundeformed.min(), self.Zundeformed.max())

        return self.Xundeformed.min(), self.Xundeformed.max(), \
               self.Yundeformed.min(), self.Yundeformed.max(), \
               self.Zundeformed.min(), self.Zundeformed.max()

    def get_limits_deformation(self):

        print("Xdeformation limits: ", self.Xdeformation.min(), self.Xdeformation.max())
        print("Ydeformation limits: ", self.Ydeformation.min(), self.Ydeformation.max())
        print("Zdeformation limits: ", self.Zdeformation.min(), self.Zdeformation.max())

        return self.Xdeformation.min(), self.Xdeformation.max(), \
               self.Ydeformation.min(), self.Ydeformation.max(), \
               self.Zdeformation.min(), self.Zdeformation.max()

    def get_limits_deformed(self):

        print("X deformed limits: ", self.Xdeformed().min(), self.Xdeformed().max())
        print("Y deformed limits: ", self.Ydeformed().min(), self.Ydeformed().max())
        print("Z deformed limits: ", self.Zdeformed().min(), self.Zdeformed().max())

        return self.Xdeformed().min(), self.Xdeformed().max(), \
               self.Ydeformed().min(), self.Ydeformed().max(), \
               self.Zdeformed().min(), self.Zdeformed().max()

    def triangulate(self):
        # triangulation
        self.triPi = numpy.array([self.Xdeformed(), self.Ydeformed()]).transpose()
        self.tri = spatial.Delaunay(self.triPi)

    def plot_triangulation(self):
        plt.triplot(self.Xdeformed(), self.Ydeformed() , self.tri.simplices.copy())
        plt.plot(self.Xdeformed(), self.Ydeformed(), "or", label = "Data")
        plt.grid()
        plt.legend()
        plt.title("triangulation")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.show()

    def get_Xinterpolated_mesh(self):
        return numpy.outer(self.x_interpolated,numpy.ones_like(self.y_interpolated))

    def get_Yinterpolated_mesh(self):
        return numpy.outer(numpy.ones_like(self.x_interpolated),self.y_interpolated)


    def interpolate(self,nx,ny,remove_nan=False):
        if self.tri is None:
            self.triangulate()

        lim = self.get_limits_deformed()
        self.x_interpolated = numpy.linspace(lim[0],lim[1],nx)
        self.y_interpolated = numpy.linspace(lim[2],lim[3],ny)

        X_INTERPOLATED =  self.get_Xinterpolated_mesh()
        Y_INTERPOLATED =  self.get_Yinterpolated_mesh()

        self.P = numpy.array([X_INTERPOLATED.flatten(), Y_INTERPOLATED.flatten() ]).transpose()

        if remove_nan:
            self.Z_INTERPOLATED = interpolate.griddata(self.triPi, self.Zdeformed(), self.P, method = "cubic", fill_value=self.Zdeformed().min() ).reshape([nx,ny])
        else:
            self.Z_INTERPOLATED = interpolate.griddata(self.triPi, self.Zdeformed(), self.P, method="cubic").reshape([nx, ny])


    def plot_interpolated(self):
        plt.contourf(self.get_Xinterpolated_mesh(), self.get_Yinterpolated_mesh(), self.Z_INTERPOLATED, 50, cmap = mpl.cm.jet)
        plt.colorbar()
        plt.contour(self.get_Xinterpolated_mesh(), self.get_Yinterpolated_mesh(), self.Z_INTERPOLATED, 20, colors = "k")
        plt.plot(self.Xdeformed(), self.Ydeformed(), "or", label = "Data")
        plt.legend()
        # plt.title = "Interpolated"  <---- THIS MAKES ERROR IN THE NEXT PLOT!!!!!!!!!!!!!!!!!
        plt.grid()
        plt.show()



    def plot_surface_image(self,invert_axes_names=True):
        if invert_axes_names:
            plot_image(self.Z_INTERPOLATED,self.x_interpolated,self.y_interpolated ,title="file: %s, axes names INVERTED from ANSYS"%self.filename,
                       xtitle="Y (%d pixels, max:%f)"%(self.x_interpolated.size,self.x_interpolated.max()),
                       ytitle="X (%d pixels, max:%f)"%(self.y_interpolated.size,self.y_interpolated.max()) )
        else:
            plot_image(self.Z_INTERPOLATED,self.x_interpolated,self.y_interpolated,title="file: %s, axes as in ANSYS"%self.filename,
                       xtitle="X (%d pixels, max:%f)"%(self.x_interpolated.size,self.x_interpolated.max()),
                       ytitle="Y (%d pixels, max:%f)"%(self.y_interpolated.size,self.y_interpolated.max()) )


    def does_interpolated_have_nan(self):
        return numpy.isnan(self.Z_INTERPOLATED).sum() > 0

    def remove_borders_in_interpolated_data(self):
        self.x_interpolated = self.x_interpolated[1:-2].copy()
        self.y_interpolated = self.y_interpolated[1:-2].copy()
        self.Z_INTERPOLATED = self.Z_INTERPOLATED[1:-2,1:-2].copy()

    def detrend(self,axis=0):
        if axis == 0:
            xm = self.x_interpolated.copy()
            zm = self.Z_INTERPOLATED[:,self.y_interpolated.size//2]
        elif axis == 1:
            xm = self.y_interpolated.copy()
            zm = self.Z_INTERPOLATED[self.x_interpolated.size // 2, :]

        zm.shape = -1

        icut = numpy.argwhere( xm > -xm[-1])
        xcut = xm[icut]
        zmcut = zm[icut]

        xcut.shape = -1
        zmcut.shape = -1

        # print(">>>>>>>>>>>>>>>>>>",xm.shape,zm.shape,xcut.shape,zmcut.shape)
        # plot(xm,zm,xcut,zmcut,title="sssss")

        print( numpy.argwhere(numpy.isnan(self.Z_INTERPOLATED)) )
        coeff = numpy.polyfit(xcut, zmcut, deg=2)

        zfit = coeff[1] * xm  + coeff[0]

        if axis ==0:
            for i in range(self.Z_INTERPOLATED.shape[1]):
                self.Z_INTERPOLATED[:,i] -= zfit
        elif axis == 1:
            for i in range(self.Z_INTERPOLATED.shape[0]):
                self.Z_INTERPOLATED[i,:] -= zfit

    def reset_height_to_minimum(self):
        self.Z_INTERPOLATED -= self.Z_INTERPOLATED.min()

    def reset_height_to_central_value(self):
        self.Z_INTERPOLATED -= self.Z_INTERPOLATED[self.Z_INTERPOLATED.shape[0]//2,self.Z_INTERPOLATED.shape[1]//2]


    def write_h5_surface(self,filename='presurface.hdf5',invert_axes_names=False):
        if invert_axes_names:
            write_generic_h5_surface(self.Z_INTERPOLATED.T, self.y_interpolated, self.x_interpolated, filename=filename)
        else:
            write_generic_h5_surface(self.Z_INTERPOLATED,self.x_interpolated,self.y_interpolated,filename=filename)


if __name__ == "__main__":
    from srxraylib.plot.gol import plot_image, set_qt, plot
    import os

    set_qt()


    o1 = FEA_File.process_file("s4.txt", n_axis_0=301, n_axis_1=51,
                 filename_out="/home/manuel/Oasys/s4.h5", invert_axes_names=True,
                 detrend=True, reset_height_method=1, do_plot=False)

    o1.plot_triangulation()
    o1.plot_interpolated()
    o1.plot_surface_image()