import h5py
from srxraylib.plot.gol import plot, plot_image, set_qt

set_qt()

def from_hrd5_to_ascii(filename_in,filename_out,do_plot=True):

    f = h5py.File(filename_in,'r')
    x = f["entry/data0/x"][:]
    y = f["entry/data0/y"][:]
    image = f["entry/data0/image"][:].T
    f.close()
    print(x.shape,y.shape,image.shape)

    if do_plot:
        plot_image(image,x,y,aspect='auto',xtitle="X / mm",ytitle="Y / mm",title="Power density W/mm2")

    print("Total power",image.sum()*(x[1]-x[0])*(y[1]-y[0]))

    f = open(filename_out,'w')
    for i in range(image.shape[0]):
        for j in range(image.shape[1]):
            f.write("%g %g %g\n"%(1e-3*x[i],1e-3*y[j],image[i,j]))

    f.close()
    print("File written to disk: %s"%filename_out)


# from_hrd5_to_ascii("C:/Users/Manuel/Oasys/cosmic_M1_H.h5","cosmic_M1_H.txt")

def extract_from_nexus(filename_in):


    f = h5py.File(filename_in,'r')

    x = f["entry/data0/x"][:]
    y = f["entry/data0/y"][:]
    z = f["entry/data0/image"][:].T

    f.close()

    return z,x,y


z,x,y = extract_from_nexus("C:/Users/Manuel/Oasys/cosmic_M1_H_xoppy.h5")
za,xa,ya = extract_from_nexus("C:/Users/Manuel/Oasys/cosmic_M1_H_arriving_xoppy.h5")

nx = x.size
ny = y.size

plot_image(z,x,y,aspect='auto')

plot(y,z[nx//2,:],
     ya,za[nx//2,:],
     ya,30.78/54.962*za[nx//2,:],
     legend=["absorbed","incident","incident times averg absorp"],
     xtitle="Y / mm",ytitle="Power density W/mm2")

plot(x,z[:,ny//2],
     xa,za[:,ny//2],
     xa,30.78/54.962*za[:,ny//2],
     legend=["absorbed","incident","incident times averg absorp"],
     xtitle="X / mm",ytitle="Power density W/mm2")