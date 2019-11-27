import numpy
from srxraylib.plot.gol import plot


def eq6(u):
    return numpy.sqrt(u) * (3 * u - 1) - 3 * p0 / 2 / numpy.abs(bm_radius) * (1 - u**2) * (1 - u)


if __name__ == "__main__":


    U = numpy.linspace(0.2,6.0,1000)


    # Moreno
    p0 = 12 * U / (1 + U)
    bm_radius = 5.28
    # ALS-U
    # p0 = 5.87 * U / (1 + U)
    # bm_radius = 7.6136

    F = eq6(U)

    F1 = F.copy()
    F1[(F1.size // 3):-1] = 10000
    imin = (numpy.abs(F1)).argmin()
    print("Zero at uvalue: %f, M= %f"%(U[imin],1.0/U[imin]))


    FF = F.copy()
    FF[0:FF.size//3] = 1000
    imin = (numpy.abs(FF)).argmin()
    print("Zero at uvalue: %f, M= %f" % (U[imin], 1.0 / U[imin]))


    plot(U, F, yrange=[-5, 5])