import numpy


def func_ellipse(x, p, q, theta):
    """
    calculates the ellipse y(x) defined by its distance to focii (p,q) and grazing
    angle theta at coordinate x=0
    :param x: the length coordinate for the ellipse (x=0 is the center)
    :param p: the distance from source to mirror center
    :param q: the distance from mirror center to image
    :param theta: the grazing incidence angle in rad
    :return:
    """


    a = (p + q) / 2

    b = numpy.sqrt( numpy.abs(p * q)) * numpy.sin(theta)

    c = numpy.sqrt(numpy.abs(a*a - b*b))

    epsilon = c / a

    # (x0,y0) are the coordinates of the center of the mirror
    # x0 = (p*p - q*q) / 4 / c
    x0 = (p - q) / 2 / epsilon
    y0 = -b * numpy.sqrt(numpy.abs(1.0 - ((x0/a)**2)))

    # the versor normal to the surface at the mirror center is -grad(ellipse)
    xnor = -2 * x0 / a**2
    ynor = -2 * y0 / b**2
    modnor = numpy.sqrt(xnor**2 + ynor**2)
    xnor /= modnor
    ynor /= modnor
    # tangent  versor is perpendicular to normal versor
    xtan =  ynor
    ytan = -xnor

    A = 1/b**2
    B = 1/a**2
    C = A

    CCC = numpy.zeros(11) # this is for the general 3D case (we need 10 coeffs, index=0 not used here)
    # The 2D implicit ellipse equation is c2 x^2 + c3 y^2 + c5 x y + c8 x + c9 y + c10 = 0
    #CCC[1] = A
    CCC[2] = B*xtan**2 + C*ytan**2
    CCC[3] = B*xnor**2 + C*ynor**2
    #CCC[4] = .0
    CCC[5] = 2*(B*xnor*xtan+C*ynor*ytan)
    #CCC[6] = .0
    #CCC[7] = .0
    CCC[8] = .0
    CCC[9] = 2*(B*x0*xnor+C*y0*ynor)
    CCC[10]= .0

    #reorder in y and get the second degree equation for heights
    # AA y^2 + BB y + CC = 0
    AA = CCC[3]
    BB = CCC[5]*x + CCC[9]
    CC = CCC[2]*x*x + CCC[8]*x + CCC[10]
    DD = BB*BB-4*AA*CC
    y = (-BB - numpy.sqrt(DD)) / 2 / AA


    return y


def ken_equation(x, p, q, theta):

    s = numpy.sin(theta)
    c = numpy.cos(theta)

    y = s * (p + q) / ( c**2 * (p - q)**2 + 4 * p * q)
    y *= 2 * p * q - c * (p - q) * x - 2 * numpy.sqrt(p * q) * numpy.sqrt(p * q - c * (p - q) * x - x**2)

    return y

def func_ellipse_heights_amparo(x1, p, q, theta, hshift=0.0, vshift=0.0):
    #
    # returns y(x), the heights of an ellipse defined by p,q, and theta using the formula in http://doi.org/10.1117/12.736171
    #

    x = x1 + hshift

    a = (p + q) / 2
    b = numpy.sqrt( p * q) * numpy.sin(theta)
    c = numpy.sqrt(a*a - b*b)

    epsilon = c / a

    # (x0,y0) are the coordinates of the center of the mirror
    # x0 = (p*p - q*q) / 4 / c
    x0 = (p - q) / 2 / epsilon
    # y0 = (-b * numpy.sqrt(1 - ((x0/a)**2)))
    # mu = numpy.arctan( - b * x0 / numpy.sqrt(a*a - x0*x0) )
    #from Amparo Excel
    # y0 = -0.99994
    # mu = 0.003894
    y0 = -b*numpy.sqrt(1-x0*x0/a/a)
    alpha = numpy.arcsin(p/2/c*numpy.sin(numpy.pi-2*theta))
    mu = alpha - theta

    brk0 = ( (x*numpy.cos(mu)+x0)/a )**2
    brk1 = b*b * (1 - brk0 )
    brk = -numpy.sqrt(brk1) -y0
    pnc = numpy.cos(mu) * ( y0 + b*numpy.sqrt( 1 - x0*x0/a/a) )
    yell2 = numpy.cos(mu) * brk + numpy.sin(-mu) * x * numpy.cos(mu) + pnc
    return yell2+vshift


def func_ellipse_heights_xianbo(x1, S1, S2, TH, hshift=0.0, vshift=0.0): # OFFSET, XC):
    #
    # Hi Manuel,
    # I have checked the mirror 20 and 21. I got quite good elliptical fit with my Igor fitting function. The figure error results are:
    # Mirror 20: sigma_s = 0.42 urad, sigma_h = 2.06 nm
    # Mirror 21: sigma_s = 0.50 urad, sigma_h = 6.62 nm
    #
    # The 1-D elliptical fitting function I used is:
    #
    # Variables: OFFSET, S1, S2 and TH, XC (x offset)
    #
    # y = OFFSET+ cos(TH-asin((S1*sin(2.0*TH))/sqrt(S1^2.0+S2^2.0+2*S1*S2*cos(2*TH))))*(sqrt(S1*S2*(1.0-(-2.0*S2*sqrt((S2+S1*cos(2*TH))^2.0/(S1^2.0+S2^2.0+2.0*S1*S2*cos(2.0*TH)))+sqrt(S1^2.0+S2^2.0+2.0*S1*S2*cos(2.0*TH)))^2.0/(S1+S2)^2.0)*sin(TH)^2.0)-sqrt(S1*S2*(1.0-(1.0/(S1+S2)^2.0)*(-2.0*S2*sqrt((S2+S1*cos(2.0*TH))^2.0/(S1^2.0+S2^2.0+2.0*S1*S2*cos(2.0*TH)))+sqrt(S1^2.0+S2^2.0+2.0*S1*S2*cos(2.0*TH))+2.0*(X-XC)*cos(TH-asin((S1*sin(2.0*TH))/sqrt(S1^2.0+S2^2.0+2.0*S1*S2*cos(2.0*TH)))))^2.0)*sin(TH)^2.0)+(X-XC)*sin(TH-asin((S1*sin(2.0*TH))/sqrt(S1^2.0+S2^2.0+2.0*S1*S2*cos(2*TH)))))

    X = x1 + hshift
    XC = 0.0
    y =         numpy.cos(TH-numpy.arcsin((S1*numpy.sin(2.0*TH))/numpy.sqrt(S1**2.0+S2**2.0+2*S1*S2*numpy.cos(2*TH))))*\
                (numpy.sqrt(S1*S2*(1.0-(-2.0*S2*numpy.sqrt((S2+S1*numpy.cos(2*TH))**2.0/ \
                (S1**2.0+S2**2.0+2.0*S1*S2*numpy.cos(2.0*TH)))+numpy.sqrt(S1**2.0+S2**2.0+ \
                2.0*S1*S2*numpy.cos(2.0*TH)))**2.0/(S1+S2)**2.0)*numpy.sin(TH)**2.0)- \
                numpy.sqrt(S1*S2*(1.0-(1.0/(S1+S2)**2.0)*(-2.0*S2*numpy.sqrt((S2+S1*numpy.cos(2.0*TH))**2.0/ \
                (S1**2.0+S2**2.0+2.0*S1*S2*numpy.cos(2.0*TH)))+numpy.sqrt(S1**2.0+S2**2.0+ \
                2.0*S1*S2*numpy.cos(2.0*TH))+2.0*(X-XC)*numpy.cos(TH-numpy.arcsin((S1*numpy.sin(2.0*TH))/ \
                numpy.sqrt(S1**2.0+S2**2.0+2.0*S1*S2*numpy.cos(2.0*TH)))))**2.0)*numpy.sin(TH)**2.0)+ \
                (X-XC)*numpy.sin(TH-numpy.arcsin((S1*numpy.sin(2.0*TH))/numpy.sqrt(S1**2.0+S2**2.0+2.0*S1*S2*numpy.cos(2*TH)))))

    return y+vshift


if __name__ == "__main__":
    p = 10.0
    q = 3.0
    theta = 3e-3
    x = 2.4

    print("result eq Ken              ", ken_equation(x, p, q, theta))
    print("result eq SHADOW+Manolo    ", func_ellipse(x, p, q, theta))
    print("result eq Amparo Rommeveaux", func_ellipse_heights_amparo(x, p, q, theta))
    print("result eq Xianbo           ", func_ellipse_heights_xianbo(x, p, q, theta))

    # result eq Ken               0.003073887585011368
    # result eq SHADOW+Manolo     0.003073887585011367
    # result eq Amparo Rommeveaux 0.003073906715650554
    # result eq Xianbo            0.0030739067156503787


    # x = numpy.linspace(-1.0, 1.0, 100)
    # from srxraylib.plot.gol import plot, set_qt
    # set_qt()
    # plot(x, ken_equation(x, p, q, theta),
    #      x, func_ellipse(x, p, q, theta), legend=["ken", "manolo"])
    #
    # print("Max diff: ", (numpy.abs(ken_equation(x, p, q, theta) - func_ellipse(x, p, q, theta))).max() )
