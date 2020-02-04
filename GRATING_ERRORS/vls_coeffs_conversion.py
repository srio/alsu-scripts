import numpy
from numpy.polynomial.polynomial import Polynomial
from srxraylib.plot.gol import plot




def from_howells_to_density(d0=0.0, v1=0.0, v2=0.0, v3=0.0, v4=0.0):
    # https://www.wolframalpha.com/input/?i=Series%281%2F%28+d0+*+%281%2Bv1+x+%2Bv2+x%5E2+%2B+v3+x%5E3%29%29%29
    s0 = 1 / d0
    s1 = - v1 / d0
    s2 = (v1 **2 - v2) / d0
    s3 = -(v1**3 - 2 * v1 * v2 + v3) / d0
    s4 = ( v1 **4 - 3 * v1**2 * v2 + 2 * v1 * v3 + v2 **2) / d0

    print("d0: %g, s0: %g" % (d0, s0))
    print("v1: %g, s1: %g" % (v1, s1))
    print("v2: %g, s2: %g" % (v2, s2))
    print("v3: %g, s3: %g" % (v3, s3))
    print("v4: %g, s4: %g" % (v4, s4))

    return s0,s1,s2,s3,s4

def from_density_to_howells(s0=0.0,s1=0.0,s2=0.0,s3=0.0,s4=0.0):

    # https://www.wolframalpha.com/input/?i=Series%281%2F%28++s0%2Bs1+x+%2Bs2+x%5E2+%2B+s3+x%5E3%29%29
    d0 = 1/s0
    v1 = - s1 / s0**2 / d0
    v2 = (s1**2 - s0 * s2) / s0**3 / d0
    v3 = -(s0 ** 2 * s3 - 2 * s0 * s1 * s2 + s1 **3) / s0**4 / d0
    v4 = ( 2 * s0**2 * s1 * s3 + s0**2 * s2**2 - 3 * s0 * s1**2 * s2 + s1**4) / s0**5 / d0

    print("d0: %g, s0: %g" % (d0, s0))
    print("v1: %g, s1: %g" % (v1, s1))
    print("v2: %g, s2: %g" % (v2, s2))
    print("v3: %g, s3: %g" % (v3, s3))
    print("v4: %g, s4: %g" % (v4, s4))

    return d0,v1,v2,v3, v4


def n_coefficients_howells(d0,v1,v2=0,v3=0, v4=0, v5=0):

    n100 = 1/d0
    n200 = -v1 / 2 / d0
    n300 = (v1**2 - v2) / 3 / d0
    n400 = (-v1**3 + 2 * v1 * v2 - v3) / 4 / d0
    n500 = (v1**4 +-3 * v1**2 * v2 + v2**2 + 2 * v1 * v3 - v4) / 5 / d0
    n600 = (-v1**5 + 4 * v1**3 * v2 - 3 * v1**2 * v3 + 2 * v2 * v3 - 3 * v1 * v2**2 + 2 * v1 * v4 - v5) / 6 / d0

    return [n100,n200,n300,n400,n500,n600]


if __name__ == "__main__":


    # print("from_howells_to_density")
    # s0,s1,s2,s3,s4 = from_howells_to_density(0.000166667,-0.001333333,4.444444e-7,0,0)
    # print("\n\nfrom_density_to_howells")
    # d0, v1, v2, v3, v4 = from_density_to_howells(s0,s1,s2,s3,s4)


    theta_inc_deg = 88.794647
    theta_ref_deg = -87.912865
    theta_grazing_inc = (90 - theta_inc_deg) * numpy.pi / 180
    theta_grazing_ref = (90 - theta_ref_deg) * numpy.pi / 180

    order = +1
    s0 = 287440.0
    s1 = 125705.19002201
    s2 = 21808.27740956
    s3 = 4004.40547258
    s4 = 0.0

    d0, v1, v2, v3, v4 = from_density_to_howells(s0,s1,s2,s3,s4)

    x = numpy.linspace(-0.003,0.003,100)
    xx = x / numpy.sin(theta_grazing_inc)

    vtot = d0 * ( 1 + v1 * xx + v2 * xx ** 2 ) #+ v3 * xx ** 3)
    stot = s0 + s1 * xx + s2 * xx**2 #+ s3 * xx**3 + s4 * xx**4
    #
    plot(xx, stot, xx, 1/vtot, xtitle="x [m]", ytitle="line density",legend=["density coeffs", "Howells"])

    print(n_coefficients_howells(d0, v1, v2, v3, v4))




    #
    #
    #

    # wf = in_object_1.duplicate()
    # wavelength = wf.get_wavelength()
    # print("wavelength = ", wavelength)
    # x = wf.get_abscissas()
    # k = wf.get_wavenumber()
    # print("wavenumber = ", k)
    #
    # #
    # # mirror
    # #
    # # f = 1.0 / (1.0/(15.098+9.203) + 1.0/7.573)
    # # print("f = ",f)
    # # phi = -k * x * x / 2 / f
    # # from srxraylib.plot.gol import plot
    # # plot(x,phi)
    # # wf.add_phase_shift(phi)
    #
    # #
    # # grating
    # #
    # order = +1
    # a0 = 287440.0
    # a1 = 125705.19002201
    # a2 = 21808.27740956
    # a3 = 4004.40547258
    #
    # # a2 = 0
    # # a3 = 0
    #
    #
    # theta_inc_deg = 88.794647
    # theta_ref_deg = -87.912865
    #
    # theta_grazing_inc = (90 - theta_inc_deg) * numpy.pi / 180
    # theta_grazing_ref = (90 - theta_ref_deg) * numpy.pi / 180
    #
    # theta_inc = numpy.pi / 2 - theta_grazing_inc
    # theta_ref = numpy.pi / 2 - theta_grazing_ref
    #
    # xx = x / numpy.sin(theta_grazing_inc)
    # xxx = xx * numpy.sin(theta_grazing_ref)
    #
    # atot = a0 * xx ** 0 + a1 * xx ** 1 + 0 * a2 * xx ** 2 + a3 * xx ** 3
    # aintegral = a0 * xx ** 1 + a1 / 2 * xx ** 2 + a2 / 3 * xx ** 3 + a3 / 4 * xx ** 4
    # aintegral_corrected = a0 * xx ** 1 + a1 / 2 * xx ** 2 + (
    #     a2) / 3 * xx ** 3  # +  (-a1**3 / a0**2 + 2 * a1 * a2 - a3)/4 * xx**4
    # # dspacing = 1.0 / atot
    # # d0 = 1.0 / a0
    #
    # path_diff1 = -xx * (numpy.sin(theta_inc) + numpy.sin(theta_ref))
    # path_diff2 = -xx * order * wavelength * a0
    # # plot(xxx, path_diff1,xxx, path_diff2, legend=["-x * (sin(inc) + sin(ref))","- x n lambda / d"],show=1)
    #
    # phi1 = k * path_diff1
    # phi2 = -xx * order * a0 * 2 * numpy.pi  # k * path_diff2
    #
    # # phi22 = -xx * atot  * 2 * numpy.pi * order
    # phi22 = - aintegral * 2 * numpy.pi * order
    #
    # path_in = numpy.sqrt(xx ** 2 - x ** 2)
    # path_out = numpy.sqrt(xx ** 2 - xxx ** 2)
    # path_tot = path_in - path_out
    # phi22 = path_tot * k - aintegral_corrected * 2 * numpy.pi * order
    # # axis correction
    # # phi22 += k * xx * (numpy.sin(theta_inc) + numpy.sin(theta_ref))
    #
    # # plot(xxx,phi22-phi2,title="phi diff")
    #
    #
    # # ca = wf.get_complex_amplitude()
    # # ca *= numpy.exp(1j * phi22)
    # # wf.set_complax_amplitude(ca)
    # wf.add_phase_shifts(phi22)
    #
    # print(">>>focal : ", 1 / (wavelength * order * a1), 1 / (1 / (15.095 + 9.203) + 1 / (5.573 + 2)))
    #
    # from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D
    #
    # wf2 = GenericWavefront1D.initialize_wavefront_from_range(xxx.min(), xxx.max(), xxx.size, wavelength)
    # wf2.set_complex_amplitude(numpy.flip(wf.get_complex_amplitude()))
    # out_object = wf2