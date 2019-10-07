#
# script to run the wiggler preprocessor (created by ShadowOui:Wiggler)
#
from srxraylib.sources import srfunc
from srxraylib.plot.gol import plot, plot_image
from srxraylib.sources.srfunc import sync_g1, sync_f
import numpy
import scipy.constants as codata


ener_gev=2.0
current = 0.5

(traj, pars) = srfunc.wiggler_trajectory(
    b_from=0,
    inData="",
    nPer=137,
    nTrajPoints=501,
    ener_gev=ener_gev,
    per=0.0288,
    kValue=3.07,
    trajFile="tmp.traj",
    shift_x_flag=0,
    shift_x_value=0.0,
    shift_betax_flag=0,
    shift_betax_value=0.0)

#  0      1    2      3      4     5       6         7
#L x[m]  y[m]  z[m]  BetaX  BetaY  BetaZ  Curvature  B[T]
print(traj.shape)

xp = traj[3,:]
b0 = traj[7,:]

# plot(1e6*xp,numpy.abs(b0),xtitle="xp")

ii = numpy.argsort(xp)

XP = numpy.linspace(xp.min()*1.2,xp.max()*1.2,400)
B0 = numpy.interp(XP,xp[ii],numpy.abs(b0[ii]))

plot(1e6*XP,B0,xtitle="XP")
# plot(xp,traj[7,:],xtitle="x'",ytitle="B")

# critical energy
Ec = 665.0 * ener_gev**2 * numpy.abs(B0)

print("Ec from %f eV to %f eV"%(Ec.min(), Ec.max()))
plot(1e6*XP,Ec,xtitle="x' / urad",ytitle="Ec / eV")






#
# flux
#
calculate_spectrum = True

enerMin = 100.0 # 100.0,
enerMax = 15000.0 # 10000.0,
enerN = 500

ENERGIES = numpy.linspace(enerMin,enerMax,enerN)

e, f, w = srfunc.wiggler_spectrum(traj,
                                  enerMin=enerMin,
                                  enerMax=enerMax,
                                  nPoints=enerN,
                                  electronCurrent=current,
                                  outFile="spectrum.dat",
                                  elliptical=False)

if False:

    FLUX = numpy.zeros_like(ENERGIES)
    POWER = numpy.zeros_like(ENERGIES)
    for i,energy in enumerate(ENERGIES):
        #TODO: chech this funny numbers....
        flux_vs_theta = 27.45 * 2.457e17 * ener_gev * current * sync_g1(energy/Ec)
        flux =  numpy.trapz(flux_vs_theta,XP) # flux_vs_theta.sum() * (XP[1]-XP[0]) #
        FLUX[i] = flux
        POWER[i] = flux * 1e3 * codata.e
        print(e[i],w[i],energy,POWER[i],w[i]/POWER[i])

    plot(ENERGIES,POWER,e,w,xtitle="Energy / eV",ytitle="W/eV",legend=["integrated","reference"])




#
# 2D theta energy
#




if False:

    plot(1e6 * XP, sync_g1(100.0 / Ec),
         1e6 * XP, sync_g1(1000.0 / Ec),
         1e6 * XP, sync_g1(10000.0 / Ec),
         xtitle="x' / eV", ytitle="Intensity", legend=["100 eV", "1000 eV ", "10000 eV"])

    THETAENERGY = numpy.zeros((XP.size, ENERGIES.size))

    for i,energy in enumerate(ENERGIES):
        #TODO: chech this funny numbers....
        flux_vs_theta = 27.45 * 2.457e17 * ener_gev * current * sync_g1(energy/Ec)
        power_vs_theta = flux_vs_theta * 1e3 * codata.e
        THETAENERGY[:, i] = power_vs_theta


    plot_image(THETAENERGY.T,ENERGIES,1e6*XP,aspect='auto',xtitle="Energy / eV",ytitle="x' / urad",show=False)
    plot(1e6 * XP, THETAENERGY.sum(axis=1) * (ENERGIES[1] - ENERGIES[0]), xtitle="x'", ytitle="Intensity", show=False)
    plot(ENERGIES, THETAENERGY.sum(axis=0) * (XP[1] - XP[0]),
         e,w,xtitle="Energy / eV",legend=["integrated","reference"], show=True)

    # plot(1e6 * XP, THETAENERGY[:,-1],
    #      1e6 * XP, 27.45 * 2.457e17 * ener_gev * current * sync_g1(ENERGIES[-1]/Ec),
    #      xtitle="x'", ytitle="Intensity latest energy %f eV"%ENERGIES[-1], show=True)


#
# psi energy
#

psi_max = 800e-6 #10 / gamma
PSI = numpy.linspace(-psi_max,psi_max,100)
codata_mee = 1e-6 * codata.m_e * codata.c**2 / codata.e
gamma = ener_gev*1e3/codata_mee


if False:
    PSITHETA = numpy.zeros((PSI.size,XP.size))


    POWER = numpy.zeros_like(ENERGIES)
    for i,energy in enumerate(ENERGIES):
        PSITHETA = sync_f(PSI * gamma, rEnergy=energy / Ec, polarization=0, gauss=0, l2=1, l3=0)
        PSITHETA /= PSITHETA.sum() * (PSI[1] - PSI[0]) * (XP[1] - XP[0])
        PSITHETA *= w[i]
        POWER[i] = PSITHETA.sum() * (PSI[1] - PSI[0]) * (XP[1] - XP[0])


    plot_image(PSITHETA.T,1e6*XP,1e6*PSI,aspect='auto',title="Energy: %f eV "%energy,ytitle="Psi / urad",xtitle="theta / urad")
    plot(ENERGIES,POWER) #,e,w,legend=["integrated","reference"])





#
# 3D energy,theta,psi
#

if True:
    ENERGYPSITHETA = numpy.zeros( (ENERGIES.size,PSI.size,XP.size))

    print("ENERGYPSITHETA.shape: ",ENERGYPSITHETA.shape)
    for i,energy in enumerate(ENERGIES):

        PSITHETA = sync_f(PSI * gamma, rEnergy=energy / Ec, polarization=0, gauss=0, l2=1, l3=0)
        PSITHETA /= PSITHETA.sum() * (PSI[1] - PSI[0]) * (XP[1] - XP[0])
        PSITHETA *= w[i]

        ENERGYPSITHETA[i,:,:] = PSITHETA



    plot_image(ENERGYPSITHETA.sum(axis=0),1e6*PSI,1e6*XP,aspect='auto',ytitle="XP",xtitle="PSI",show=False)

    plot(1e6*XP,ENERGYPSITHETA.sum(axis=0).sum(axis=0),xtitle="X'",show=False)

    plot(1e6*PSI,ENERGYPSITHETA.sum(axis=0).sum(axis=1),xtitle="PSI")

