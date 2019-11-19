from srxraylib.plot.gol import plot_image
import h5py
import matplotlib.pylab as plt
import numpy

def run_xoppy():
    #
    # script to make the calculations (created by XOPPY:undulator_spectrum)
    #
    from orangecontrib.xoppy.util.xoppy_undulators import xoppy_calc_undulator_power_density

    h5_parameters = dict()
    h5_parameters["ELECTRONENERGY"] = 2.0
    h5_parameters["ELECTRONENERGYSPREAD"] = 0.00104
    h5_parameters["ELECTRONCURRENT"] = 0.5
    h5_parameters["ELECTRONBEAMSIZEH"] = 1.212e-05
    h5_parameters["ELECTRONBEAMSIZEV"] = 1.473e-05
    h5_parameters["ELECTRONBEAMDIVERGENCEH"] = 5.77e-06
    h5_parameters["ELECTRONBEAMDIVERGENCEV"] = 4.75e-06
    h5_parameters["PERIODID"] = 0.0288
    h5_parameters["NPERIODS"] = 137.0
    h5_parameters["KV"] = 3.07
    h5_parameters["KH"] = 0.0
    h5_parameters["KPHASE"] = 0.0
    h5_parameters["DISTANCE"] = 13.73
    h5_parameters["GAPH"] = 0.03
    h5_parameters["GAPV"] = 0.015
    h5_parameters["HSLITPOINTS"] = 201
    h5_parameters["VSLITPOINTS"] = 201
    h5_parameters["METHOD"] = 2
    h5_parameters["USEEMITTANCES"] = 1
    h5_parameters["MASK_FLAG"] = 0
    h5_parameters["MASK_ROT_H_DEG"] = 0.0
    h5_parameters["MASK_ROT_V_DEG"] = 88.75
    h5_parameters["MASK_H_MIN"] = -50.0
    h5_parameters["MASK_H_MAX"] = 50.0
    h5_parameters["MASK_V_MIN"] = -1.2
    h5_parameters["MASK_V_MAX"] = 1.2

    h, v, p, code = xoppy_calc_undulator_power_density(
        ELECTRONENERGY=h5_parameters["ELECTRONENERGY"],
        ELECTRONENERGYSPREAD=h5_parameters["ELECTRONENERGYSPREAD"],
        ELECTRONCURRENT=h5_parameters["ELECTRONCURRENT"],
        ELECTRONBEAMSIZEH=h5_parameters["ELECTRONBEAMSIZEH"],
        ELECTRONBEAMSIZEV=h5_parameters["ELECTRONBEAMSIZEV"],
        ELECTRONBEAMDIVERGENCEH=h5_parameters["ELECTRONBEAMDIVERGENCEH"],
        ELECTRONBEAMDIVERGENCEV=h5_parameters["ELECTRONBEAMDIVERGENCEV"],
        PERIODID=h5_parameters["PERIODID"],
        NPERIODS=h5_parameters["NPERIODS"],
        KV=h5_parameters["KV"],
        KH=h5_parameters["KH"],
        KPHASE=h5_parameters["KPHASE"],
        DISTANCE=h5_parameters["DISTANCE"],
        GAPH=h5_parameters["GAPH"],
        GAPV=h5_parameters["GAPV"],
        HSLITPOINTS=h5_parameters["HSLITPOINTS"],
        VSLITPOINTS=h5_parameters["VSLITPOINTS"],
        METHOD=h5_parameters["METHOD"],
        USEEMITTANCES=h5_parameters["USEEMITTANCES"],
        MASK_FLAG=h5_parameters["MASK_FLAG"],
        MASK_ROT_H_DEG=h5_parameters["MASK_ROT_H_DEG"],
        MASK_ROT_V_DEG=h5_parameters["MASK_ROT_V_DEG"],
        MASK_H_MIN=h5_parameters["MASK_H_MIN"],
        MASK_H_MAX=h5_parameters["MASK_H_MAX"],
        MASK_V_MIN=h5_parameters["MASK_V_MIN"],
        MASK_V_MAX=h5_parameters["MASK_V_MAX"],
        h5_file="undulator_power_density.h5",
        h5_entry_name="XOPPY_POWERDENSITY",
        h5_initialize=True,
        h5_parameters=h5_parameters,
    )


    return p,h,v,

def get_from_h5(filename):

    f = h5py.File(filename,'r')

    p = f["XOPPY_POWERDENSITY/PowerDensity/image_data"][:].T
    h = f["XOPPY_POWERDENSITY/PowerDensity/axis_x"][:]
    v = f["XOPPY_POWERDENSITY/PowerDensity/axis_y"][:]

    f.close()

    return p,h,v



def integrate_power_density(p,h,v):
    hplus = numpy.array(numpy.where(h >= 0))
    hplus.shape = -1
    vplus = numpy.array(numpy.where(v >= 0))
    vplus.shape = -1
    p2 = p[hplus[0]:1 + hplus[-1], vplus[0]:1 + vplus[-1]].copy()
    pintegrated = numpy.cumsum(p2, 0)
    pintegrated = numpy.cumsum(pintegrated, 1)
    pintegrated *= numpy.abs((h[1] - h[0]) * (v[1] - v[0])) * 4  # the four is because we double the interval in H and V

    hintegrated = 2.0 * h[hplus]  # axes are now gap (aperture)
    vintegrated = 2.0 * v[vplus]  # axes are now gap (aperture)

    return pintegrated, hintegrated, vintegrated




def plot_image_contour(p,h,v,aspect='equal',title="",xtitle="",ytitle="",show=1):

    fig = plt.figure()

    # cmap = plt.cm.Greys
    plt.imshow(p.T,origin='lower',extent=[h[0],h[-1],v[0],v[-1]],cmap=None,aspect=aspect)
    if True:
        plt.colorbar()
    ax = fig.gca()
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)

    plt.title(title)

    levels = numpy.arange(0.0, p.max()*0.95, 100.0)  # Boost the upper limit to avoid truncation errors.

    vv = plt.axis()
    ff = plt.contour(p.T, levels, colors='k', origin='lower', extent=[h[0],h[-1],v[0],v[-1]])
    plt.clabel(ff, fmt='%d', colors='b', fontsize=14)
    plt.axis(vv)

    if show:
        plt.show()



if __name__ == "__main__":
    # example plot

    # p, h, v = run_xoppy()
    # plot_image(p, h, v, xtitle="H [mm]", ytitle="V [mm]", title="Power density W/mm2")

    p, h, v = get_from_h5("/home/manuel/Oasys/undulator_power_density.h5")
    plot_image(p, h, v, xtitle="H [mm]", ytitle="V [mm]", title="Power density W/mm2", show=0)

    pintegrated, hintegrated, vintegrated = integrate_power_density(p,h,v)



    nh = int(0.5*hintegrated.size)
    nv = int(0.5*vintegrated.size)
    fig3 = plot_image_contour(pintegrated[0:nh,0:nv],hintegrated[0:nh],vintegrated[0:nv],aspect='equal',title="Power through slit [W]",
                              xtitle="Slit gap X[mm]",ytitle="Slit gap Y [mm]",show=1)

