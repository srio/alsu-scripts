from srxraylib.plot.gol import plot
import matplotlib.pylab as plt
import numpy


import matplotlib
matplotlib.rc('xtick',         labelsize=20)
matplotlib.rc('ytick',         labelsize=20)
matplotlib.rcParams.update({'font.size': 20})


# ap = numpy.loadtxt( "flexon_ken_memo2_factor1_fit.dat", skiprows=1)

# is_fit = True
#
# if is_fit:
#     is_fit_string = "_fit"
# else:
#     is_fit_string = ""

photon_energy1 = 250
photon_energy2 = 1250

for is_fit_string in ["_fit",""]:
    am = numpy.loadtxt("flexon_ken_memo2_cosine%s.dat"%(is_fit_string), skiprows=1)
    fm =    plot(am[:,0], am[:,1],
                 am[:,0], am[:,2],
                 am[:,0], am[:,3],
                 am[:,0], am[:,4],
                 xlog=False,figsize=(12,8),
                 legend=["Uncorrected E=%d eV"%photon_energy1,"Corrected E=%d eV"%photon_energy1,"Uncorrected E=%d eV"%photon_energy2,"Corrected E=%d eV"%photon_energy2],
                 xtitle="Number of ripples in mirror length",ytitle="Strehl I/I0",
                 show=0)

    matplotlib.pylab.grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.2)
    # matplotlib.pylab.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)

    filename = "flexon_ken_memo2_cosine%s.png"%(is_fit_string)
    fm[0].savefig(filename)

    matplotlib.pylab.show()
    print("File %s written to disk"%filename)


