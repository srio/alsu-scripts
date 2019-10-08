import numpy
import scipy.constants as codata

import xraylib
from srxraylib.plot.gol import plot, plot_scatter
m2ev = codata.c * codata.h / codata.e      # lambda(m)  = m2eV / energy(eV)

cosd = lambda x : numpy.cos( numpy.deg2rad(x) )
sind = lambda x : numpy.sin( numpy.deg2rad(x) )
asind = lambda x : numpy.rad2deg( numpy.arcsin( x ))
acosd = lambda x : numpy.rad2deg( numpy.arccos( x ))

# # used by Luca
# codata_h = numpy.array(6.62606957e-34)
# codata_ec = numpy.array(1.602176565e-19)
# codata_c = numpy.array(299792458.0)
# m2ev = codata_c * codata_h / codata_ec      # lambda(m)  = m2eV / energy(eV)

def solve_grating_equation(line_density=100000,wavelength=10e-10,c=1.10,order=1,method='beta'):

    if method == 'beta':  # beta!!

        A = (1.0 - 1.0 / c ** 2)
        B = -2 * order * line_density * wavelength
        C = order ** 2 * line_density ** 2 * wavelength ** 2 - 1 + 1.0 / c ** 2

        Delta = B ** 2 - 4 * A * C

        sinbeta1 = (-B + numpy.sqrt(Delta)) / 2 / A
        sinbeta2 = (-B - numpy.sqrt(Delta)) / 2 / A

        # print("Discriminant=%f, sinbeta1=%f, sinbeta2=%f" % (Delta, sinbeta1, sinbeta2))

        if numpy.abs(sinbeta1) <= 1:
            sinbeta = sinbeta1
        elif numpy.abs(sinbeta2) <= 1:
            sinbeta = sinbeta2
        else:
            raise Exception("No valid beta angle")

        # change sign of beta
        # sinbeta *= -1

        # beta = numpy.arcsin(sinbeta)

        sinalpha = order * wavelength * line_density - sinbeta

        # alpha = numpy.arcsin(sinalpha)

    else:  # alpha
        A = (1.0 - order ** 2)
        B = -2 * order * line_density * wavelength
        C = order ** 2 * line_density ** 2 * wavelength ** 2 - 1 + c ** 2

        Delta = B ** 2 - 4 * A * C

        sinalpha1 = (-B + numpy.sqrt(Delta)) / 2 / A
        sinalpha2 = (-B - numpy.sqrt(Delta)) / 2 / A

        print("Discriminant=%f, sinalpha1=%f, sinalpha2=%f" % (Delta, sinalpha1, sinalpha2))

        if numpy.abs(sinalpha1) <= 1:
            sinalpha = sinalpha1
        elif numpy.abs(sinalpha2) <= 1:
            sinalpha = sinalpha2
        else:
            raise Exception("No valid alpha angle")

        # # my value
        # sinalpha_me = m*lambda0*k0/(1-c**2) + numpy.sqrt( (m*lambda0*k0/(1-c**2))**2 - \
        #                     ( (m * lambda0 * k0 / (1-c**2))**2 -1) )
        #
        # # Luca
        # sinalpha_luca = (-m * k0 * lambda0 / (c ** 2 - 1)) + \
        #             numpy.sqrt(1 + (m * m * c * c * k0 * k0 * lambda0 * lambda0) / (
        #                         (c ** 2 - 1) ** 2))
        #
        # print("sin_alpha numeric, me, luca: ",sinalpha,sinalpha_me,sinalpha_luca)

        # change sign of beta
        # sinbeta *= -1

        # alpha = numpy.arcsin(sinalpha)

        sinbeta = order * wavelength * line_density - sinalpha

        # beta = numpy.arcsin(sinbeta)

    return sinalpha,sinbeta

def interface_reflectivity(alpha,gamma,theta1):
    """
    Calculates the reflectivity of an interface using Fresnel formulas.

    Code adapted from XOP and SHADOW
    :param alpha: the array with alpha values (alpha=2 delta, n=1-delta+i beta)
    :param gamma: the array with gamma (gamma=2 beta)
    :param theta1: a scalar with the grazing angle in rad
    :return:
    """

    rho =  numpy.sin(theta1)**2 - alpha
    rho += numpy.sqrt( ( (numpy.sin(theta1))**2 - alpha)**2 + gamma**2)
    rho = numpy.sqrt(rho / 2)
    # ;** Computes now the reflectivity for s-pol


    rs1 = 4 * (rho**2) * (numpy.sin(theta1) - rho)**2 + gamma**2
    rs2 = 4 * (rho**2) * (numpy.sin(theta1) + rho)**2 + gamma**2
    rs = rs1 / rs2

    # ;** Computes now the polarization ratio


    ratio1 = 4 * rho**2 * (rho * numpy.sin(theta1) - numpy.cos(theta1))**2 + gamma**2 * numpy.sin(theta1)**2
    ratio2 = 4 * rho**2 * (rho * numpy.sin(theta1) + numpy.cos(theta1))**2 + gamma**2 * numpy.sin(theta1)**2
    ratio = ratio1 / ratio2

    rp = rs * ratio
    runp = 0.5 * (rs + rp)

    return rs,rp,runp

def reflectivity(descriptor,energy,theta,density=None,rough=0.0):



    if isinstance(descriptor,str):
        Z = xraylib.SymbolToAtomicNumber(descriptor)
        symbol = descriptor
    else:
        Z = descriptor
        symbol = xraylib.AtomicNumberToSymbol(descriptor)

    if density == None:
        density = xraylib.ElementDensity(Z)

    atwt = xraylib.AtomicWeight(Z)
    avogadro = codata.Avogadro
    toangstroms = codata.h * codata.c / codata.e * 1e10
    re = codata.e**2 / codata.m_e / codata.c**2 / (4*numpy.pi*codata.epsilon_0) * 1e2 # in cm

    molecules_per_cc = density * avogadro / atwt
    wavelength = toangstroms / energy  * 1e-8 # in cm
    k = molecules_per_cc * re * wavelength * wavelength / 2.0 / numpy.pi

    f1 = numpy.zeros_like(energy)
    f2 = numpy.zeros_like(energy)
    for i,ienergy in enumerate(energy):
        f1[i] = Z + xraylib.Fi(Z,1e-3*ienergy)
        f2[i] = - xraylib.Fii(Z,1e-3*ienergy)

    alpha = 2.0 * k * f1
    gamma = 2.0 * k * f2

    rs,rp,runp = interface_reflectivity(alpha,gamma,theta)

    if rough != 0:
        rough *= 1e-8 # to cm
        debyewaller = numpy.exp( -( 4.0 * numpy.pi * numpy.sin(theta) * rough / wavelength)**2)
    else:
        debyewaller = 1.0

    return rs*debyewaller,rp*debyewaller,runp*debyewaller

def structure_factor(lamda,gamma_blaze,alpha,beta,d):
    gamma_blaze_deg = gamma_blaze * 180.0 / numpy.pi
    alphadash = 90 - alpha * 180 / numpy.pi
    betadash = 90 + beta * 180 / numpy.pi
    b = d * (numpy.sin(numpy.pi / 180 * alphadash) / numpy.sin(numpy.pi / 180 * (alphadash + gamma_blaze_deg)))
    factor1 = numpy.pi * b / lamda
    factor2 = numpy.cos(numpy.pi / 180 * (alphadash + gamma_blaze_deg)) - numpy.cos(numpy.pi / 180 * (betadash - gamma_blaze_deg))
    factor3 = factor1 * factor2
    s = (numpy.sin(factor3)) / factor3
    return s

def run_jark_2019():
    # inputs
    energy = numpy.linspace(250,4000,300)
    g_inv_mm = 150.0
    C = 1.245
    # gamma_blaze = 0.20 * numpy.pi / 180
    energy_blaze = 800.0

    #
    wavelength = codata.h * codata.c / codata.e / energy
    # grating
    g_inv_m = g_inv_mm * 1e3
    dspacing = (1.0 / g_inv_m)

    print("dspacing: %f A "%(1e10*dspacing))

    alpha = numpy.zeros_like(wavelength)
    beta = numpy.zeros_like(wavelength)

    for i in range(wavelength.size):
        sinalpha, sinbeta = solve_grating_equation(line_density=g_inv_m, wavelength=wavelength[i], c=C, order=1, method='beta')
        alpha[i] = numpy.arcsin(sinalpha)
        beta[i] =  numpy.arcsin(sinbeta)

    theta = (alpha - beta) / 2


    # plot(energy,alpha*180/numpy.pi,
    #      energy, -beta * 180 / numpy.pi,
    #      energy, (alpha-beta) / 2 * 180 / numpy.pi,
    #      xtitle="energy / eV", ytitle="deg", legend=["alpha","-beta","theta"])

    RS = numpy.zeros_like(energy)

    theta = 0.5 * (alpha - beta)
    gamma = 0.5 * (alpha + beta)

    gamma_blaze = numpy.interp(energy_blaze,energy,gamma)
    print(">>>> using blaze angle: %f deg "%(gamma_blaze*180/numpy.pi))

    #
    # reflectivity
    #
    for i in range(energy.size):
        rs, rp, runp = reflectivity("Au", numpy.array([energy[i]]), numpy.pi / 2 - theta[i], )
        RS[i] = rs[0]
    # plot(energy, RS, xlog=True, xtitle="energy / eV",ytitle="R(theta)")

    #
    # efficiency
    #

    efficiency = RS / C

    sf = structure_factor(wavelength,gamma_blaze,alpha,beta,dspacing)

    efficiency *= sf**2

    # plot(energy, sf**2, xlog=True, xtitle="energy / eV", ytitle="S")

    plot(energy,efficiency,xlog=True,title='C = %f , N = %d l/mm '%(C,g_inv_mm),xtitle="Phonon energy [eV]",ytitle="Efficiency" )

def run_shadow():
    #
    # Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
    #
    import Shadow
    import numpy

    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 0

    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = Shadow.Beam()
    oe0 = Shadow.Source()
    oe1 = Shadow.OE()
    oe2 = Shadow.OE()
    oe3 = Shadow.OE()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.FDISTR = 3
    oe0.F_COLOR = 3
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.0
    oe0.HDIV2 = 0.0
    oe0.ISTAR1 = 5676561
    oe0.NPOINT = 50000
    oe0.PH1 = 805.5
    oe0.PH2 = 806.5
    oe0.SIGDIX = 1.9733083328114146e-05
    oe0.SIGDIZ = 1.946778306932498e-05
    oe0.SIGMAX = 1.7218556115924652e-05
    oe0.SIGMAZ = 1.913527305050143e-05
    oe0.VDIV1 = 0.0
    oe0.VDIV2 = 0.0

    oe1.ALPHA = 90.0
    oe1.DUMMY = 100.0
    oe1.FWRITE = 1
    oe1.T_IMAGE = 0.0
    oe1.T_INCIDENCE = 88.31811450000001
    oe1.T_REFLECTION = 88.31811450000001
    oe1.T_SOURCE = 24.301

    oe2.ALPHA = 180.0
    oe2.DUMMY = 100.0
    oe2.FWRITE = 1
    oe2.F_ANGLE = 1
    oe2.F_GRATING = 1
    oe2.F_RULING = 5
    oe2.F_RUL_ABS = 1
    oe2.RULING = 300000.0
    oe2.RUL_A1 = 131195.47246228
    oe2.RUL_A2 = 22760.09661057
    oe2.RUL_A3 = 4178.99967999
    oe2.T_IMAGE = 7.573
    oe2.T_INCIDENCE = 88.768556
    oe2.T_REFLECTION = 87.867673
    oe2.T_SOURCE = 0.0

    oe3.ALPHA = 90.0
    oe3.DUMMY = 100.0
    oe3.FCYL = 1
    oe3.FMIRR = 2
    oe3.FWRITE = 1
    oe3.F_DEFAULT = 0
    oe3.SIMAG = 5.573
    oe3.SSOUR = 26.301
    oe3.THETA = 88.75
    oe3.T_IMAGE = 5.573
    oe3.T_INCIDENCE = 88.75
    oe3.T_REFLECTION = 88.75
    oe3.T_SOURCE = -5.573

    # Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")

    #
    # run optical element 1
    #
    print("    Running optical element: %d" % (1))
    if iwrite:
        oe1.write("start.01")

    beam.traceOE(oe1, 1)

    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")

    #
    # run optical element 2
    #
    print("    Running optical element: %d" % (2))
    if iwrite:
        oe2.write("start.02")

    beam.traceOE(oe2, 2)

    if iwrite:
        oe2.write("end.02")
        beam.write("star.02")

    #
    # run optical element 3
    #
    print("    Running optical element: %d" % (3))
    if iwrite:
        oe3.write("start.03")

    beam.traceOE(oe3, 3)

    if iwrite:
        oe3.write("end.03")
        beam.write("star.03")

    # Shadow.ShadowTools.plotxy(beam, 3, 1, nbins=101, nolost=1, title="Real space")

    return beam


def calculate_efficiency(mirr=None,angle_file="angle.02",material="Au",gamma_blaze_deg=None,lines_per_mm=300.0):

    import Shadow

    beam = Shadow.Beam()
    beam.load(mirr)

    angle = numpy.loadtxt(angle_file)

    alpha_deg = angle[:, 1]
    beta_deg = -angle[:, 2]
    wavelength = beam.getshonecol(19) * 1e-10
    energy = beam.getshonecol(11)

    # plot(numpy.arange(alpha_deg.size),alpha_deg,
    #      numpy.arange(alpha_deg.size), -beta_deg,
    #      legend=["alpha","beta"])
    # print(">>>",angle.shape)

    # plot(energy,wavelength)

    theta_deg = (alpha_deg - beta_deg) / 2
    gamma_deg = (alpha_deg + beta_deg) / 2

    # plot(energy,alpha*180/numpy.pi,
    #      energy, -beta * 180 / numpy.pi,
    #      energy, (alpha-beta) / 2 * 180 / numpy.pi,
    #      xtitle="energy / eV", ytitle="deg", legend=["alpha","-beta","theta"])

    RS = numpy.zeros_like(energy)

    if gamma_blaze_deg is None:
        gamma_blaze_deg = gamma_deg.mean()
        print(">>>> using blaze angle: %f deg " % (gamma_blaze_deg))

    #
    # reflectivity
    #
    theta = theta_deg * numpy.pi / 180

    method = 1

    if method == 0:
        for i in range(energy.size):
            rs, rp, runp = reflectivity("Au", numpy.array([energy[i]]), numpy.pi / 2 - theta[i], )
            RS[i] = rs[0]
    else:
        n2_interp = numpy.zeros(energy.size, complex)
        for i in range(energy.size):
            n2_interp[i] = xraylib.Refractive_Index(material, 1e-3 * energy[i],
                                                    xraylib.ElementDensity(xraylib.SymbolToAtomicNumber(material)))
        n1 = 1 + 0j  # % air
        k0 = 2 * numpy.pi / wavelength
        # % calculate the k values
        k1 = (n2_interp ** 2 - 1 + (sind(90.0 - theta_deg)) ** 2) ** 0.5
        # % calculate fresnel coefficients for the 0-1 and 1-2 interfaces
        r1 = (sind(90.0 - theta_deg) - k1) / (sind(90.0 - theta_deg) + k1)
        RS = (numpy.abs(r1)) ** 2

    plot_scatter(energy, RS, xtitle="energy / eV", ytitle="R(theta)")

    #
    # efficiency
    #

    C = numpy.cos(numpy.deg2rad(beta_deg)) / numpy.cos(numpy.deg2rad(alpha_deg))

    plot_scatter(energy, C)
    efficiency = RS / C

    dspacing = 1e-3 / lines_per_mm
    sf = structure_factor(wavelength,
                          numpy.deg2rad(gamma_blaze_deg),
                          numpy.deg2rad(alpha_deg),
                          numpy.deg2rad(beta_deg),
                          dspacing)

    efficiency *= sf ** 2

    plot_scatter(energy, sf ** 2, xtitle="energy / eV", ytitle="S**2", title="S**2")
    plot_scatter(energy, efficiency, title="", xtitle="Phonon energy [eV]", ytitle="Efficiency")

    print("results: ")
    print("Reflectivity R min: %f, max: %f, average: %f" % (RS.min(), RS.max(), RS.mean()))
    print("C min: %f, max: %f, average: %f" % (C.min(), C.max(), C.mean()))
    print("sf**2 min: %f, max: %f, average: %f" % ((sf ** 2).min(), (sf ** 2).max(), (sf ** 2).mean()))
    print("Efficiency min: %f, max: %f, average: %f" % (efficiency.min(), efficiency.max(), efficiency.mean()))

    return efficiency

if __name__ == "__main__":

    # beam = run_shadow()

    eff = calculate_efficiency(mirr="mirr.02",angle_file="angle.02",material="Au",gamma_blaze_deg=0.45,lines_per_mm=300.0)
    plot_scatter(numpy.arange(eff.size),eff,title="Efficiency",xtitle="ray index)",ytitle="Efficiency")

