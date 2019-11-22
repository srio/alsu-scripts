from orangecontrib.srw.util.srw_objects import SRWData
import numpy
from wofry.propagator.wavefront2D.generic_wavefront import GenericWavefront2D
from wofry.propagator.polarization import Polarization

def generic_wfr_add_grating_error(wfr,file_name,do_plot=False):

    wavelength = wfr.get_wavelength()

    print(wfr.is_polarized())


    # x (mm) measured d-spacing error (nm)
    e = numpy.loadtxt(file_name,skiprows=1)
    # plot(e[:, 0], e[:, 1])

    alpha_deg = 88.7946477
    beta_deg = -87.912865

    alpha_grazing = (90 - alpha_deg) * numpy.pi / 180
    print(e.shape)
    x = e[:, 0] * 1e-3
    x *= numpy.sin(alpha_grazing)

    y = e[:, 1] * 1e-9

    phi = y / wavelength
    phi *= numpy.sin(alpha_deg * numpy.pi / 180) + numpy.sin(beta_deg * numpy.pi / 180)

    # plot(x,phi,xtitle="Diffraction direction perpendicular to axis [m]",ytitle="Phi grating [rad]")

    yy = wfr.get_coordinate_y()
    print("wavefront size: ",wfr.size())
    print("error file x: min:%f max: %f size: %d" %(x.min(),x.max(),x.size))
    print("wavefront yy: min:%f max: %f size: %d"%(yy.min(), yy.max(), yy.size))

    phi_interpolated = numpy.interp(yy,x,phi)

    if do_plot:
        plot(x,phi,yy,phi_interpolated,legend=["original","interpolated"],
             xtitle="Diffraction direction perpendicular to axis [m]",ytitle="Phi grating [rad]")

    PHI = numpy.zeros(wfr.size())
    for i in range(wfr.size()[0]):
        PHI[i,:] = phi_interpolated

    if do_plot:
        plot_image(PHI,wfr.get_coordinate_x(),wfr.get_coordinate_y(),title="PHI in rads")

    wfr2 = wfr.duplicate()

    wfr2.add_phase_shifts( PHI, polarization=Polarization.SIGMA)
    wfr2.add_phase_shifts( PHI, polarization=Polarization.PI)

    return wfr2

if __name__ == "__main__":
    from srxraylib.plot.gol import plot, plot_image, set_qt
    #
    file_name = "/home/manuel/Oasys/tmp.h5"
    wfr = GenericWavefront2D.load_h5_file(file_name,"wfr")

    phase = wfr.get_phase(polarization=Polarization.SIGMA)
    plot_image(phase,wfr.get_coordinate_x(),wfr.get_coordinate_y(),title="before",show=False)

    wfr2 = generic_wfr_add_grating_error(wfr,"LEG_300 l_mm.txt",do_plot=True)

    phase2 = wfr2.get_phase(polarization=Polarization.SIGMA)
    plot_image(phase2,wfr2.get_coordinate_x(),wfr2.get_coordinate_y(),title="after")

    plot_image(phase2-phase,wfr2.get_coordinate_x(),wfr2.get_coordinate_y(),title="difference")