import numpy
from srxraylib.plot.gol import plot, set_qt
from scipy.ndimage import gaussian_filter1d
import scipy.constants as codata
import srxraylib.sources.srfunc as srfunc

def get_magnetic_field_ALSU(do_plot=False,filename=""):
    L = 1605.0 #mm

    y = numpy.linspace(0,L, 2000)

    B = y * 0.0

    for i in range(y.size):
        if y[i] > 75 and y[i] < 575: B[i] = -0.876
        if y[i] > 650 and y[i] < 975: B[i] = 0.16
        if y[i] > 1030 and y[i] < 1530: B[i] = -0.8497

    # plot(y, B)
    B2 = gaussian_filter1d(B, 2.5)

    yy = y.copy()
    yy -= yy[y.size//2]
    yy *= 1e-3

    if do_plot:
        # plot(yy, B, yy, B2, legend=["original","smoothed"],xtitle="y / m",ytitle="B / T")
        plot(yy, B2, xtitle="y [m]", ytitle="B [T]")

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
        plot(y, B, xtitle="y [m]", ytitle="B [T]")

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
    print("Half-Divergence AB: ", 0.5 * (0.325) / (1e9 / codata.c * electron_energy_in_GeV/0.16) )
    print("Half-Divergence M2: ", 0.5 * (0.500) / (1e9 / codata.c * electron_energy_in_GeV/0.8497) )

    get_magnetic_field_ALSU(do_plot=True,filename="BM_multi.b")
    get_magnetic_field_ALS(do_plot=True, filename="BM_als.b")