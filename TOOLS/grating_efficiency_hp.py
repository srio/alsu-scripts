# %
# % code to calculate grating efficiency assuming constant C as
# % a function of photon energy
# %
# %  H. A. Padmore    9/27/19
# %

import numpy
import xraylib
import matplotlib.pylab as plt

cosd = lambda x : numpy.cos( numpy.deg2rad(x) )
sind = lambda x : numpy.sin( numpy.deg2rad(x) )
asind = lambda x : numpy.rad2deg( numpy.arcsin( x ))
acosd = lambda x : numpy.rad2deg( numpy.arccos( x ))


#
# inputs
#
emin= 250.0     # float( input('input the minimum photon energy  ') )
emax= 4000.0    # float( input('input the maximum photon energy  ') )
einc= 100.0     # float( input('input the increment in photon energy  ') )
c= 1.245        # float( input('input mono constant C ') )
g= 150.0        # float( input('input grating line density, l/mm ') )
energy_b= 800.0 # float( input('input blaze energy, eV ') )

d=1/(g*1000)
lamda_b=1e-10*12398.5/energy_b
const_a=1-1/c**2
const_b=-2*lamda_b/d
const_c=(lamda_b/d)**2-1+(1/c**2)
beta  = asind((-const_b-(const_b**2-4*const_a*const_c)**0.5)/(2*const_a))
alpha = asind((lamda_b/d)-sind(beta))
blaze=(alpha+beta)/2
print('Using suggested blaze angle is %4.2f °\n'%blaze)
gamma= blaze # float(input('input blaze angle, ° '))
energy= numpy.arange(emin,emax,einc)
len_ein=energy.size
# %
# % calculate diffraction efficiency in constant deviation mode
# %
# % calculate alpha and beta for constant C
lamda=1e-10*12398.5/energy
const_a=1-1/c**2
const_b=-2*lamda/d
const_c=(lamda/d)**2-1+(1/c**2)
beta= asind((-const_b-(const_b**2-4*const_a*const_c)**0.5)/(2*const_a))
alpha=asind((lamda/d)-sind(beta))
alphadash=90-alpha
betadash=90+beta
b=d*(sind(alphadash)/sind(alphadash+gamma))
factor1=numpy.pi*b/lamda
factor2=cosd(alphadash+gamma)-cosd(betadash-gamma)
factor3=factor1*factor2
s=(numpy.sin(factor3))/factor3
bare_eff=(alphadash/betadash)*s**2
#     %
#     % calculate reflectivity at the complement of the 1/2 included angle
#     %
theta=(alpha-beta)/2
thetadash=90-theta

n2_interp = numpy.zeros(energy.size,complex)
for i in range(energy.size):
    n2_interp[i] = xraylib.Refractive_Index("Au",1e-3*energy[i],
                                            xraylib.ElementDensity(xraylib.SymbolToAtomicNumber("Au")))

n1=1+0j  #% air
k0=2*numpy.pi/lamda
# % calculate the k values
k1=(n2_interp**2-1+(sind(thetadash))**2)**0.5
# % calculate fresnel coefficients for the 0-1 and 1-2 interfaces
r1=(sind(thetadash)-k1)/(sind(thetadash)+k1)
R=(numpy.abs(r1))**2

#     % calculate efficiency including reflectivity
eff_r=bare_eff*R

# plot results
plt.plot(energy,eff_r)
plt.xscale('log')
plt.xlabel("Phonon energy [eV]")
plt.ylabel("Efficiency")
plt.title('C = %f , N = %d l/mm '%(c,g))
plt.show()

