import numpy
from srxraylib.plot.gol import plot, set_qt
from scipy.ndimage import gaussian_filter1d
import scipy.constants as codata
import srxraylib.sources.srfunc as srfunc


def get_magnetic_field_ALSU_onlyMag7(do_plot=False,filename=""):

    drift = 75.0
    lengthBM = 500.0


    L = 2 * drift + lengthBM

    y = numpy.linspace(0,L, 2000)

    B = y * 0.0



    B0_7 = -0.876
    for i in range(y.size):

        if y[i] > drift and y[i] < drift+lengthBM: B[i] = B0_7

    # plot(y, B)
    B2 = gaussian_filter1d(B, 2.5)

    yy = y.copy()
    yy -= drift + lengthBM / 2
    yy *= 1e-3

    if do_plot:
        # plot(yy, B, yy, B2, legend=["original","smoothed"],xtitle="y / m",ytitle="B / T")
        plot(yy, B2, xtitle="y [m]", ytitle="B [T]",title=filename)

    if filename != "":
        f = open(filename, "w")
        for i in range(y.size):
            f.write("%f  %f\n" % (yy[i], B2[i]))
        f.close()
        print("File written to disk: %s"%filename)

    return yy,B2



def get_magnetic_field_ALSU_centeredMag7(do_plot=False,filename=""):

    drift = 75.0
    lengthBM = 500.0
    lengthAB = 305


    L = 4 * drift + 2 * lengthAB + lengthBM  # 1605.0 #mm
    L = 5 * drift + 2 * lengthAB + 2 * lengthBM  # 1605.0 #mm

    y = numpy.linspace(0,L, 2000)

    B = y * 0.0



    B0_7 = -0.876
    B0_AB = 0.16
    B0_8 = -0.8497
    for i in range(y.size):


        # if y[i] > drift and y[i] < drift+lengthBM: B[i] = -0.876
        # if y[i] > 2*drift+lengthBM and y[i] < 2*drift+lengthBM+lengthAB: B[i] = 0.16
        # if y[i] > 3*drift+lengthBM+lengthAB and y[i] < 3*drift+2*lengthBM+lengthAB: B[i] = -0.8497

        if y[i] > drift and y[i] < drift+lengthAB: B[i] = B0_AB
        if y[i] > 2*drift+lengthAB and y[i] < 2*drift+lengthAB+lengthBM: B[i] = B0_7
        if y[i] > 3*drift+lengthAB+lengthBM and y[i] < 3*drift+2*lengthAB+lengthBM: B[i] = B0_AB
        if y[i] > 4*drift+2*lengthAB+lengthBM and y[i] < 4*drift+2*lengthAB+2*lengthBM: B[i] = B0_8


    # plot(y, B)
    B2 = gaussian_filter1d(B, 2.5)

    yy = y.copy()
    yy -= 2 * drift + lengthAB + lengthBM / 2
    yy *= 1e-3

    if do_plot:
        # plot(yy, B, yy, B2, legend=["original","smoothed"],xtitle="y / m",ytitle="B / T")
        plot(yy, B2, xtitle="y [m]", ytitle="B [T]",title=filename)

    if filename != "":
        f = open(filename, "w")
        for i in range(y.size):
            f.write("%f  %f\n" % (yy[i], B2[i]))
        f.close()
        print("File written to disk: %s"%filename)

    return yy,B2


def get_magnetic_field_ALSU(do_plot=False,filename=""):


    L = 1605.0 #mm

    y = numpy.linspace(0,L, 2000)

    B = y * 0.0

    drift = 75.0
    lengthBM = 500.0
    lengthAB = 305

    for i in range(y.size):


        if y[i] > drift and y[i] < drift+lengthBM: B[i] = -0.876
        if y[i] > 2*drift+lengthBM and y[i] < 2*drift+lengthBM+lengthAB: B[i] = 0.16
        if y[i] > 3*drift+lengthBM+lengthAB and y[i] < 3*drift+2*lengthBM+lengthAB: B[i] = -0.8497


    # plot(y, B)
    B2 = gaussian_filter1d(B, 2.5)

    yy = y.copy()
    yy -= yy[y.size//2]
    yy *= 1e-3

    if do_plot:
        # plot(yy, B, yy, B2, legend=["original","smoothed"],xtitle="y / m",ytitle="B / T")
        plot(yy, B2, xtitle="y [m]", ytitle="B [T]",title=filename)

    if filename != "":
        f = open(filename, "w")
        for i in range(y.size):
            f.write("%f  %f\n" % (yy[i], B2[i]))
        f.close()
        print("File written to disk: %s"%filename)

    return yy,B2

def get_magnetic_field_ALS(do_plot=False,filename=""):

    L = 0.350
    Lmax = 0.6
    y = numpy.linspace(-Lmax / 2, Lmax / 2, 500)

    B = y * 0.0 - 1.27

    ybad = numpy.where(numpy.abs(y) > (L / 2))

    B[ybad] = 0

    if do_plot:
        plot(y, B, xtitle="y [m]", ytitle="B [T]",title=filename)

    if filename != "":
        f = open(filename, "w")
        for i in range(y.size):
            f.write("%f  %f\n" % (y[i], B[i]))
        f.close()
        print("File written to disk: %s"%filename)

    return y,B


if __name__ == "__main__":
    set_qt()
    electron_energy_in_GeV = 2.0

    print("Radius M1: ", 1e9 / codata.c * electron_energy_in_GeV/0.876)
    print("Radius AB: ", 1e9 / codata.c * electron_energy_in_GeV/0.16)
    print("Radius M2: ", 1e9 / codata.c * electron_energy_in_GeV/0.849)

    print("Half-Divergence M1: ", 0.5 * (0.500) / (1e9 / codata.c * electron_energy_in_GeV/0.876) )
    print("Half-Divergence AB: ", 0.5 * (0.305) / (1e9 / codata.c * electron_energy_in_GeV/0.16) )
    print("Half-Divergence M2: ", 0.5 * (0.500) / (1e9 / codata.c * electron_energy_in_GeV/0.8497) )

    get_magnetic_field_ALSU(do_plot=True,filename="BM_multi.b")
    get_magnetic_field_ALS(do_plot=True, filename="BM_als.b")
    get_magnetic_field_ALSU_centeredMag7(do_plot=True,filename="BM_multi7.b")
    get_magnetic_field_ALSU_onlyMag7(do_plot=True, filename="BM_only7.b")