

def get_stack_from_h5(filein,path="/source/Radiation stack"):
    f = h5py.File(filein, 'r')
    e = f["%s/axis0"%path][:]
    x = f["%s/axis1"%path][:]
    y = f["%s/axis2"%path][:]
    rad0 = f["%s/stack_data"%path][:]
    f.close()
    return e,x,y,rad0

def get_image_from_h5(filein,path=""):
    f = h5py.File(filein, 'r')
    d = f["%s/image_data"%path][:].T
    x = f["%s/axis_x"%path][:]
    y = f["%s/axis_y"%path][:]
    f.close()
    return x,y,d


import h5py

# filein = "/home/manuel/Oasys/undulator_radiation_flexon_grant3.h5"
#
# filein = "/home/manuel/Oasys/power3D.h5"


# e,x,y,rad0 = get_stack_from_h5(filein="/home/manuel/Oasys/power3D.h5",path="/source/Radiation stack")
# print(">>>>",e.shape,x.shape,y.shape,rad0.shape)
# e,x,y,rad1 = get_stack_from_h5(filein="/home/manuel/Oasys/power3D.h5",path="/optical_element_1/Transmittivity_stack")
# print(">>>>",e.shape,x.shape,y.shape,rad1.shape)

x0,y0,powerdensity0 = get_image_from_h5(filein="/home/manuel/Oasys/power3D.h5",path="/source/Power Density")
# print(">>>>",x.shape,y.shape,powerdensity0.shape)

x1,y1,powerdensity1 = get_image_from_h5(filein="/home/manuel/Oasys/power3D.h5",path="/optical_element_1/Absorbed Power Density")
# print(">>>>",x.shape,y.shape,powerdensity0.shape)

# step_x = x[1] - x[0]
# step_y = y[1] - y[0]
# step_e = e[1] - e[0]


from srxraylib.plot.gol import plot_image, plot



print("Total power at source: ",     powerdensity0.sum()*(x0[1]-x0[0])*(y0[1]-y0[0]))
print("Total power absorbed my M1: ",powerdensity1.sum()*(x1[1]-x1[0])*(y1[1]-y1[0]))

plot_image(powerdensity0,x0,y0,title="source",xtitle="X [um]",ytitle="Y [um]",show=False)
plot_image(powerdensity1,x1,y1,title="mirror",xtitle="X [um]",ytitle="Y [um]",show=False)


plot(x0,powerdensity0[:,151//2],
     x1,powerdensity1[:,151//2],xtitle="X [um]",ytitle="Power density [W/mm2]",legend=["source","absorbed M1"],show=0)

plot(y0,powerdensity0[151//2,:],
     y1,powerdensity1[151//2,:],xtitle="Y [um]",ytitle="Power density [W/mm2]",legend=["source","absorbed M1"],)
