import numpy
from srxraylib.plot.gol import plot

def create_grating(lines_per_m=300000.0,grating_length=0.150,points_per_period=9,
                   form_code=0,vls=[300000.0,269816.234363,87748.010405,27876.983114]):
    """
    Create a grating profile

    :param lines_per_m:
    :param grating_length:
    :param points_per_period:
    :param form_code:  =0sin, 1=cos, 2=square
    :return: (x,y) the arrays with abscissas and normalized ordinates (ampliture=1)

    """

    period = 1. / lines_per_m
    number_of_periods = grating_length / period
    print("Number of periods: ",number_of_periods)
    print("Period: %f um"%(1e6*period))

    x = numpy.linspace(-grating_length/2,grating_length/2,int(points_per_period * number_of_periods))
    if form_code == 0: # sin
        y = (numpy.sin(2 * numpy.pi * x / period) + 1) / 2
    elif form_code == 1: # cos
        y = (numpy.cos(2 * numpy.pi * x / period) + 1) / 2
    elif form_code == 2: # square
        from scipy.signal import square
        y = (square(2 * numpy.pi * x / period, duty=0.5) + 1) / 2
    elif form_code == 3: # vls
        from scipy.signal import sweep_poly
        p = numpy.poly1d([vls[3],vls[2],vls[1], vls[0]])
        y = numpy.ceil( sweep_poly(x, p ) )

    return x,y


if __name__ == "__main__":
    import matplotlib
    matplotlib.rc('axes.formatter', useoffset=False)
    lines_per_m=300000
    period = 1. / lines_per_m
    x,y = create_grating(lines_per_m=300000,grating_length=0.150,points_per_period=12,form_code=3)

    # y *= -1
    # y += 1
    y *= 10e-9

    # plot(1e3*x,y ,xrange=[-75,-74.98],marker='o')
    print("Number of points: ",x.size)

    filename = "C:/Users/manuel/Oasys/VLS_FLEXON_myscript.txt"
    if filename != "":
        f = open(filename,'w')
        for i in range(x.size):
            f.write("%g %g\n"%(x[i],y[i]))
        f.close()
        print("File written to disk: %s"%filename)

    b = numpy.loadtxt("C:/Users/Manuel/Oasys/VLS_FLEXON.txt")
    xb = b[:,0]
    yb = b[:,1]

    plot(1e3*x,y,1e3*xb,yb,legend=["mio","daniele"],xrange=[-75+148.98,-74.98+148.98],show=0)
    plot(1e3 * x, y, 1e3 * xb, yb, legend=["mio", "daniele"], xrange=[-75 , -74.98],show=0)
    plot(1e3 * x, y, 1e3 * xb, yb, legend=["mio", "daniele"], xrange=[-0.1, 0.1])

    # plot(1e6*x,y)
