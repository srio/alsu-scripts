import numpy
from srxraylib.plot.gol import plot, set_qt
from scipy.ndimage import gaussian_filter1d
import scipy.constants as codata
import srxraylib.sources.srfunc as srfunc


def calculate_flux(y,B,M=0):
    # analyse M1

    B3 = B.copy()

    if M==1:
        ibad = numpy.argwhere(y > -0.2)
        B3[ibad] = 0.0
    elif M==2:
        ibad = numpy.argwhere( numpy.abs(y) > 0.2)
        B3[ibad] = 0.0
    elif M==3:
        ibad = numpy.argwhere(y < 0.2)
        B3[ibad] = 0.0
    else:
        pass

    # from srxraylib.plot.gol import plot
    # plot(y,B3)

    # plot(y,B3)
    # B3[numpy.where(y > 575)] = 0.0
    #
    # # plot(yy,B2)
    # radius = 3.334728 * electron_energy_in_GeV / B
    # w = B3 / B3.min()
    # # plot(yy,w)
    # center = numpy.average(yy, weights=w)
    #
    # t = numpy.where(w > 0.5)
    # width = yy[t[0][-1]] - yy[t[0][0]]
    #
    # # print(">>>>Div M1: ",(575-75)/2000*1.605/radius.min())
    # # print(">>>>Half Div M1: ",0.5*(575-75)/2000*1.605/radius.min())
    # print(">>>>B min: ", B.min())
    # plot(yy, B3, title="radius")
    #
    # print("M1: center: %f m, width: %f m, half-divergence=%f rad" % (center, width, 0.5 * width / radius.min()))

    # f = open("BM_tmp.b", "w")
    # for i in range(y.size):
    #     f.write("%f  %f\n" % (yy[i], B3[i]))
    # f.close()
    # print("File written to disk: BM_tmp.b")

    tmp = numpy.vstack((y,B3)).T
    print(">>>>",tmp.shape,tmp[:,0])

    (traj, pars) = srfunc.wiggler_trajectory(
        b_from=1,
        inData=tmp, #"BM_first.b",
        nPer=1,
        nTrajPoints=501,
        ener_gev=1.9,
        per=0.01,
        kValue=1.0,
        trajFile="tmp.traj",
        shift_x_flag=5,
        shift_x_value=0.042,
        shift_betax_flag=4,
        shift_betax_value=0.0324)

    calculate_spectrum = True

    if calculate_spectrum:
        e, f, w = srfunc.wiggler_spectrum(traj,
                                          enerMin=0.0010,
                                          enerMax=10000.1,
                                          nPoints=500,
                                          electronCurrent=400.0 * 1e-3,
                                          outFile="spectrum.dat",
                                          elliptical=False)
        from srxraylib.plot.gol import plot

    return e,f,w


if __name__ == "__main__":
    set_qt()

    L = 1605.0 #mm
    electron_energy_in_GeV = 1.9


    y = numpy.linspace(0,L, 2000)

    B = y * 0.0


    for i in range(y.size):
        if y[i] > 75 and y[i] < 575: B[i] = -0.876
        if y[i] > 650 and y[i] < 975: B[i] = 0.16
        if y[i] > 1030 and y[i] < 1530: B[i] = -0.8497

    # plot(y, B)


    B2 = gaussian_filter1d(B, 2.5)

    yy = y.copy()
    yy -= yy[y.size//2]
    yy *= 1e-3

    # plot(yy, B, yy, B2, legend=["original","smoothed"],xtitle="y / m",ytitle="B / T")

    f = open("BM_multi.b", "w")
    for i in range(y.size):
        f.write("%f  %f\n" % (yy[i], B2[i]))
    f.close()
    print("File written to disk: BM_multi.b")






    print("Radius M1: ", 1e9 / codata.c * electron_energy_in_GeV/0.876)
    print("Radius AB: ", 1e9 / codata.c * electron_energy_in_GeV/0.16)
    print("Radius M2: ", 1e9 / codata.c * electron_energy_in_GeV/0.8497)

    print("Half-Divergence M1: ", 0.5 * (0.500) / (1e9 / codata.c * electron_energy_in_GeV/0.876) )
    print("Half-Divergence AB: ", 0.5 * (0.325) / (1e9 / codata.c * electron_energy_in_GeV/0.16) )
    print("Half-Divergence M2: ", 0.5 * (0.500) / (1e9 / codata.c * electron_energy_in_GeV/0.8497) )



    #
    # spectrum
    #

    e1, f1, w1 = calculate_flux(yy, B2, M=1)
    e2, f2, w2 = calculate_flux(yy, B2, M=2)
    e3, f3, w3 = calculate_flux(yy, B2, M=3)
    e0, f0, w0 = calculate_flux(yy, B2, M=0)

    plot(e0, f0,
         e1, f1,
         e2, f2,
         e3, f3,
         xlog=True, ylog=True, show=True, yrange=[1e11,1e16],
         xtitle="Photon energy [eV]", ytitle="Flux [Photons/s/0.1%bw]", title="Flux",
         legend=["M7+AB+M8","M7","AB","M8"])
    # plot(e, w, xlog=False, ylog=False, show=True,
    #      xtitle="Photon energy [eV]", ytitle="Spectral Power [E/eV]", title="Spectral Power")
    #


    # #
    # # script to run the wiggler preprocessor (created by ShadowOui:Wiggler)
    # #
    # from srxraylib.sources import srfunc
    #
    # # (traj, pars) = srfunc.wiggler_trajectory(
    # #     b_from=1,
    # #     inData="BM_multi.b",
    # #     nPer=1,
    # #     nTrajPoints=501,
    # #     ener_gev=1.9,
    # #     per=0.01,
    # #     kValue=1.0,
    # #     trajFile="tmp.traj",
    # #     shift_x_flag=4,
    # #     shift_x_value=0.0,
    # #     shift_betax_flag=4,
    # #     shift_betax_value=0.094)
    #
    # (traj, pars) = srfunc.wiggler_trajectory(
    #     b_from=1,
    #     inData="BM_multi.b",
    #     nPer=1,
    #     nTrajPoints=501,
    #     ener_gev=1.9,
    #     per=0.01,
    #     kValue=1.0,
    #     trajFile="tmp.traj",
    #     shift_x_flag=5,
    #     shift_x_value=0.042,
    #     shift_betax_flag=4,
    #     shift_betax_value=0.0324)
    #
    # #
    # # calculate cdf and write file for Shadow/Source
    # #
    #
    #
    # # analyse M1
    #
    # B3 = B2.copy()
    # B3[numpy.where(y > 575)] = 0.0
    #
    # # plot(yy,B2)
    # radius = 3.334728*electron_energy_in_GeV/B
    # w = B3 / B3.min()
    # # plot(yy,w)
    # center = numpy.average(yy,weights=w)
    #
    # t = numpy.where(w > 0.5)
    # width = yy[t[0][-1]]-yy[t[0][0]]
    #
    #
    # # print(">>>>Div M1: ",(575-75)/2000*1.605/radius.min())
    # # print(">>>>Half Div M1: ",0.5*(575-75)/2000*1.605/radius.min())
    # print(">>>>B min: ",B.min())
    # plot(yy,B3,title="radius")
    #
    # print("M1: center: %f m, width: %f m, half-divergence=%f rad"%(center,width,0.5*width/radius.min()))
    #
    #
    # f = open("BM_first.b", "w")
    # for i in range(y.size):
    #     f.write("%f  %f\n" % (yy[i], B3[i]))
    # f.close()
    # print("File written to disk: BM_first.b")
    #
    #
    # (traj, pars) = srfunc.wiggler_trajectory(
    #     b_from=1,
    #     inData="BM_first.b",
    #     nPer=1,
    #     nTrajPoints=501,
    #     ener_gev=1.9,
    #     per=0.01,
    #     kValue=1.0,
    #     trajFile="tmp.traj",
    #     shift_x_flag=5,
    #     shift_x_value=0.042,
    #     shift_betax_flag=4,
    #     shift_betax_value=0.0324)
    #
    #
    # calculate_spectrum = True
    #
    # if calculate_spectrum:
    #     e, f, w = srfunc.wiggler_spectrum(traj,
    #                                       enerMin=0.0010,
    #                                       enerMax=10000.1,
    #                                       nPoints=500,
    #                                       electronCurrent=400.0 * 1e-3,
    #                                       outFile="spectrum.dat",
    #                                       elliptical=False)
    #     from srxraylib.plot.gol import plot
    #
    #     plot(e, f, xlog=True, ylog=True, show=True,
    #          xtitle="Photon energy [eV]", ytitle="Flux [Photons/s/0.1%bw]", title="Flux")
    #     # plot(e, w, xlog=False, ylog=False, show=True,
    #     #      xtitle="Photon energy [eV]", ytitle="Spectral Power [E/eV]", title="Spectral Power")
    #     #