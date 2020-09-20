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

if __name__ == "__main__":
    p = 10.0
    q = 3.0
    theta = 3e-3
    x = 2.4

    print(ken_equation(x, p, q, theta))
    print(func_ellipse(x, p, q, theta))


    x = numpy.linspace(-1.0, 1.0, 100)
    from srxraylib.plot.gol import plot, set_qt
    set_qt()
    plot(x, ken_equation(x, p, q, theta),
         x, func_ellipse(x, p, q, theta), legend=["ken", "manolo"])
