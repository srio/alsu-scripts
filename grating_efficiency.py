import numpy
import scipy.constants as codata

import xraylib

m2ev = codata.c * codata.h / codata.e      # lambda(m)  = m2eV / energy(eV)


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

        # print("Discriminant=%f, sinalpha1=%f, sinalpha2=%f" % (Delta, sinalpha1, sinalpha2))

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

def structure_factor(alpha,beta,gamma,wavelength,dspacing):
    b = dspacing * numpy.sin(alpha) / numpy.sin(alpha + gamma)
    x = (numpy.pi / b / wavelength) * ( numpy.cos(alpha + gamma) - numpy.cos(beta - gamma) )
    S = numpy.sin( x ) / x
    return S

if __name__ == "__main__":
    from srxraylib.plot.gol import plot

    energy = numpy.linspace(200,4000,300)


    #
    wavelength = codata.h * codata.c / codata.e / energy
    # print("wavelengt: ",wavelength*1e10)


    #
    # grating
    #
    g_inv_mm = 150.0
    g_inv_m = g_inv_mm * 1e3
    dspacing = (1.0 / g_inv_m)

    print("dspacing: %f A "%(1e10*dspacing))
    # gamma = 0.22 * numpy.pi / 180.0
    # deflection_angle_deg = 3.0

    C = 1.245

    alpha = numpy.zeros_like(wavelength)
    beta = numpy.zeros_like(wavelength)

    for i in range(wavelength.size):
        sinalpha, sinbeta = solve_grating_equation(line_density=g_inv_m, wavelength=wavelength[i], c=C, order=1, method='beta')
        alpha[i] = numpy.arcsin(sinalpha)
        beta[i] =  numpy.arcsin(sinbeta)

    theta = (alpha - beta) / 2



    # grating  n lambda = 0.5 * dspacing (theta**2 - phi**2)

    # theta = 0.5 * (alpha - beta)  # half-include
    # gamma = 0.5 * (alpha + beta)  # blaze


    # alpha = numpy.sqrt( wavelength/2/dspacing * (C+1) / (C-1))

    # beta = 2 * gamma - alpha

    # beta = numpy.arcsin(wavelength / dspacing - alpha)


    plot(energy,alpha*180/numpy.pi,
         energy, -beta * 180 / numpy.pi,
         energy, (alpha-beta) / 2 * 180 / numpy.pi,
         xtitle="energy / eV", ytitle="deg", legend=["alpha","-beta","theta"])



    RS = numpy.zeros_like(energy)

    theta = 0.5 * (alpha - beta)

    for i in range(energy.size):
        rs, rp, runp = reflectivity("Au", numpy.array([energy[i]]), numpy.pi / 2 - theta[i], )
        # print(">>> rs: ",rs,rs.shape)
        RS[i] = rs[0]

    plot(energy, RS, xlog=True, xtitle="energy / eV",ytitle="R(theta)")

    efficiency = RS / C

    # plot(energy, efficiency, xlog=True, xtitle="energy / eV", ytitle="Eff")

    gamma = 0.5 * (alpha + beta)
    sf = structure_factor(numpy.pi / 2 - alpha,
                          numpy.pi / 2 - beta,
                          gamma,wavelength,dspacing)

    plot(energy, sf, xlog=True, xtitle="energy / eV", ytitle="S")

    # theta = 2 * deflection_angle_deg
    # gamma = blaze_angle_deg
    # alpha = theta + gamma




