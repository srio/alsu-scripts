import numpy
from srxraylib.plot.gol import plot

def create_cos(length=0.150,number_of_ripples=1.0,points_per_period=12,amplitude=20e-9):

    x = numpy.linspace(-length/2,length/2,int(points_per_period * number_of_ripples))

    period = length / number_of_ripples
    y = numpy.cos(2 * numpy.pi * x / period)

    y -= y.min()
    y /= y.max()
    y *= amplitude

    return x,y


if __name__ == "__main__":
    import matplotlib
    matplotlib.rc('axes.formatter', useoffset=False)

    x,y = create_cos(length=0.150,number_of_ripples=10.0,points_per_period=120)


    print("Number of points: ",x.size)

    filename = "C:/Users/manuel/Oasys/cosine.txt"
    if filename != "":
        f = open(filename,'w')
        for i in range(x.size):
            f.write("%g %g\n"%(x[i],y[i]))
        f.close()
        print("File written to disk: %s"%filename)


    # plot(1e3*x,y,show=1)

