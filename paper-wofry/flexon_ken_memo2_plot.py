from srxraylib.plot.gol import plot
import matplotlib.pylab as plt
import numpy


import matplotlib
matplotlib.rc('xtick',         labelsize=20)
matplotlib.rc('ytick',         labelsize=20)
matplotlib.rcParams.update({'font.size': 20})


ap = numpy.loadtxt( "flexon_ken_memo2_factor1.dat", skiprows=1)


photon_energy1 = 250
photon_energy2 = 1250

for factor in ["-1","1"]:
    am = numpy.loadtxt("flexon_ken_memo2_factor%s.dat"%factor, skiprows=1)
    fm =    plot(am[:,0], am[:,1],
                 am[:,0], am[:,2],
                 am[:,0], am[:,3],
                 am[:,0], am[:,4],
                 xlog=True,figsize=(12,8),
                 legend=["Uncorrected E=%d eV"%photon_energy1,"Corrected E=%d eV"%photon_energy1,"Uncorrected E=%d eV"%photon_energy2,"Corrected E=%d eV"%photon_energy2],
                 xtitle="Radius [m]",ytitle="Strehl I/I0",xrange=[50,1e6])

    filename = "flexon_ken_memo2_factor%s.png"%factor
    fm[0].savefig(filename)
    print("File %s written to disk"%filename)


