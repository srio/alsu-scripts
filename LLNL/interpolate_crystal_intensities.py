

import numpy
import Shadow

import h5py
from srxraylib.plot.gol import plot_image, plot, plot_scatter
from scipy import interpolate



#
# xop file
#

filename = "/home/manuel/Oasys/Ge111_spol.h5"

f = h5py.File(filename,'r')
ang=f["/ang"][:]
en=f["/en"][:]
ref=f["/ref"][:]
f.close()

# print(en.shape,ang.shape,ref.shape)
# plot_image(ref,en,ang,aspect='auto',xtitle="energy [eV]",ytitle="angle [urad]")


#
# Shadow files
#
beam = Shadow.Beam()
beam.load("/home/manuel/Oasys/star.01")
aa = numpy.loadtxt("/home/manuel/Oasys/angle.01")
energies = beam.getshonecol(11)
wavelength_cm = 2 * numpy.pi / beam.rays[:,10]
print(aa.shape)
incident_angles_deg = 90.0 - aa[:,1]
indices = aa[:,0]
Shadow.ShadowTools.plotxy(beam,11,3)



#
# calculate bragg angle for each ray energy
#


d_spacing_cm = 3.266272e-8
theta = numpy.arcsin( wavelength_cm / 2 / d_spacing_cm )

theta_deg = theta * 180 / numpy.pi

print("ratio", wavelength_cm / 2 / d_spacing_cm)
print("theta_deg: ",theta_deg)
print("wavelength_cm: ",wavelength_cm)
print("incident_angles_deg",incident_angles_deg)
print("theta_deg",theta_deg)

plot_scatter(incident_angles_deg,theta_deg,xtitle="incident angle [deg]",ytitle="Bragg angle [deg]")

print(energies.shape,incident_angles_deg.shape)
# plot(energies,incident_angles_deg)

angle_diff_in_urad = (incident_angles_deg - theta_deg) * numpy.pi / 180 * 1e6
plot_scatter(indices,angle_diff_in_urad,xtitle="index",ytitle="angle diff")


# f = interpolate.interp2d(en, ang, ref.T, kind='cubic')
# f(energies,angle_diff_in_urad)


interpolated_weight = numpy.exp( - angle_diff_in_urad**2 / 2/ (400/2.35)**2)
plot_scatter(angle_diff_in_urad,interpolated_weight)

beam.rays[:, 6] =  numpy.sqrt(interpolated_weight) * beam.rays[:, 6]
beam.rays[:, 7] =  numpy.sqrt(interpolated_weight) * beam.rays[:, 7]
beam.rays[:, 8] =  numpy.sqrt(interpolated_weight) * beam.rays[:, 8]
beam.rays[:, 15] = 0  # * beam.rays[:, 15]
beam.rays[:, 16] = 0  # * beam.rays[:, 16]
beam.rays[:, 17] = 0  # * beam.rays[:, 17]


Shadow.ShadowTools.plotxy(beam,11,3,ref=23)
