import numpy
from srxraylib.plot.gol import plot, plot_image, set_qt
from srxraylib.metrology.dabam import dabam
import h5py

from run_wofry import run_beamline, run_beamline_2

def extract_profile(filename, index, filename_out="", renormalize_to_heights_sd=None):
    file = h5py.File(filename, 'r')
    subgroup_name = "profiles_stack"
    image_x = file[subgroup_name + "/X"][()]
    x0 = file[subgroup_name + "/Y"][()]
    PROFILES = file[subgroup_name + "/Z"][()].T
    print(PROFILES.shape, image_x.shape, x0.shape)
    y0 = PROFILES[index, :].copy()

    if renormalize_to_heights_sd is not None:
        y0 = y0 / y0.std() * renormalize_to_heights_sd

    if filename_out != "":
        f = open(filename_out, 'w')
        for i in range(x0.size):
            f.write("%g %g \n" % (x0[i], y0[i]))
        f.close()
        print("Filee written to disk: %s" % filename_out)
    return x0, y0, image_x[index]

def fit_beta(filename="tmp.dat"):
    dm = dabam.initialize_from_external_data(filename,
                                             column_index_abscissas=0,
                                             column_index_ordinates=1,
                                             skiprows=0,
                                             useHeightsOrSlopes=0,
                                             to_SI_abscissas=1.0,
                                             to_SI_ordinates=1.0,
                                             detrending_flag=-1)
    # dm.plot("heights")
    # print(">>>Dabam: ",dm.y.size)
    # plot(dm.y, dm.zHeights)
    # print(dm.info_profiles())
    return dm, -dm.powerlaw["hgt_pendent"]

if __name__ == "__main__":

    set_qt()

    profile_index = 0

    # x, y, beta = extract_profile("set1_profiles.h5", profile_index, filename_out="tmp.dat",
    #                              renormalize_to_heights_sd=1e-9)
    #
    # print(">>>>>>>>>>>",y.size,y.std(),)
    #
    # dm, betafit = fit_beta()
    # zz = numpy.loadtxt("tmp.dat")
    #
    # plot(#x, y,
    #      dm.y,dm.zHeights,
    #      zz[:,0], zz[:,1], show=1, title="beta: %f" % (betafit))




    rms = numpy.linspace(0,20,21)
    out = numpy.zeros((10,rms.size))

    for i in range(rms.size):

        # filename_out = "tmp.dat"
        # x, y, beta = extract_profile("set1_profiles.h5", profile_index, filename_out=filename_out,
        #                              renormalize_to_heights_sd=rms[i] * 1e-9)


        filename_out = "C:/Users/Manuel/OASYS1.2/alsu-scripts/STREHL_TW/fractalprofile.dat"
        a = numpy.loadtxt(filename_out)
        x = a[:, 0]
        y = a[:, 1]
        print(">>>>>>>>>", i, y.std(), rms[i] * 1e-9)
        y = y / y.std() * rms[i] * 1e-9
        # plot(x,y)
        beta = 1.0
        f = open("tmp.dat", 'w')
        for ii in range(y.size):
            f.write("%g  %g\n" % (x[ii], y[ii]))
        f.close()
        print("File written to disk: tmp.dat")


        wf = run_beamline(error_flag=1,error_file="tmp.dat")
        wf2 = run_beamline_2(error_flag=1, error_file="tmp.dat")

        out[0, i] = rms[i]
        out[1, i] = wf.get_intensity().max()
        out[2, i] = wf2.get_intensity().max()

        if rms.size == 1:
            plot(1e6 * wf.get_abscissas(), wf.get_intensity())
            plot(x,  1e9 * y,
                 xtitle="x [m]", ytitle="height [nm]", title="profile #%d" % profile_index)

    delta_phi = 2 * rms * 1e-9 * numpy.sin(1.25 * numpy.pi / 180)
    sr2 = numpy.exp(-(2 * numpy.pi / wf.get_wavelength() * delta_phi) ** 2)

    dm, betafit = fit_beta(filename_out)


    plot(out[0,:], out[1,:] / out[1,0] ,
         out[0, :], out[2, :] / out[2, 0],
         rms, sr2,
         legend=["numeric (on-axis propagator)","numeric (grazing-propagator)","analytical"],
         title="Beta: %4.2f fit: %4.2f" % (beta, betafit),
         xtitle="Height Error rms [nm]", ytitle="Strehl Ratio",
         yrange=[0,1.1], show=0)
    import matplotlib
    matplotlib.pyplot.grid()
    matplotlib.pyplot.show()
    print(out)