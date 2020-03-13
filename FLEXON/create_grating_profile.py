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
        from scipy.signal import square
        period0 = 1.0 / vls[0] * numpy.ones_like(x)
        period = 1.0 / (vls[0] + vls[1] * x + vls[2] * x**2 + vls[3] * x**3)
        y = square(2 * numpy.pi * x / period, duty=(0.5 * period / period0))
        # y = square(2 * numpy.pi * x / period, duty=0.5)

    return x,y


if __name__ == "__main__":
    import matplotlib
    matplotlib.rc('axes.formatter', useoffset=False)
    lines_per_m=300000
    period = 1. / lines_per_m
    x,y = create_grating(lines_per_m=300000,grating_length=0.150,points_per_period=11,form_code=3)

    plot(1e3*x,y ,xrange=[-75,-74.98],marker='o')
    print("Number of points: ",x.size)
    # plot(1e6*x,y)
