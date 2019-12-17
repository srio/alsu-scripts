import numpy
from srxraylib.plot.gol import plot, set_qt
from scipy.optimize import fsolve

def eq6(u,totd1,bm_radius1):
    p0 = totd1 * u / (1 + u)
    return numpy.sqrt(u) * (3 * u - 1) - 3 * p0 / 2 / numpy.abs(bm_radius1) * (1 - u**2) * (1 - u)

def eq4(u,totd1,bm_radius1):
    p0 = totd1 * u / (1 + u)
    return 180 / numpy.pi * numpy.arctan( (1.0 - u**2) * 3 * p0 / 2 / bm_radius1)


if __name__ == "__main__":


    set_qt()

    for setup in range(6):

        # setup = 0 # 0=Moreno, 1=ALS, 2=ALSU Mag7, 3=ALSU-Mag8

        if setup == 0:
            # Moreno
            name="Moreno"
            totd = 12.0
            bm_radius = 5.28
            uvalue_in_use = 1.9
        elif setup == 1:
            # ALS
            name="ALS"
            totd = 5.87
            bm_radius = 5.0
            uvalue_in_use = 0.37
        elif setup == 2:
            name="ALSU-Mag7"
            # ALSU Mag7
            totd = 5.87
            bm_radius = 7.6736
            uvalue_in_use = 0.37
        elif setup == 3:
            # ALSU Mag8
            name="ALSU-Mag8"
            totd = 5.87
            bm_radius = 7.858
            uvalue_in_use = 0.37
        elif setup == 4:
            # Moreno
            name="MorenoOtherSolution"
            totd = 12.0
            bm_radius = 5.28
            uvalue_in_use = 0.52
        elif setup == 5:
            # ALSU Optimum toroid
            name="ToroidOptimum"
            totd = 5.87
            bm_radius = 7.6736
            uvalue_in_use = 3.007984

        tmp = fsolve(eq6, args=(totd,bm_radius), x0=(0.4,2.7) ) #, args=(1.2, 6))
        print("=========================== %s \n"%name)
        print("Solutions for setup %d: u1: %f, u2: %f"%(setup,tmp[0],tmp[1]))
        print("                        M1: %f, M2: %f" % (1./tmp[0], 1./tmp[1]))
        print("Angle for u in use (u=%f) grazing: %f deg, from normal: %f deg" % (uvalue_in_use,
                                                      eq4(uvalue_in_use, totd, bm_radius),
                                                      90-numpy.abs(eq4(uvalue_in_use, totd, bm_radius))
                                                      ))
        print("p: %f,  q: %f"%(totd*uvalue_in_use/(1+uvalue_in_use),totd*(1-uvalue_in_use/(1+uvalue_in_use))))

        U = numpy.linspace(0.2,6.0,1000)
        F = eq6(U,totd,bm_radius)
        plot(U, F, yrange=[-5, 5],title="setup %s: u1: %f, u2: %f"%(name,tmp[0],tmp[1]))


