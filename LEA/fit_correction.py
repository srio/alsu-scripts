import numpy



def fit_correction(filein, fileout="", calculate=False):
    from axo import orthonormalize_a, linear_2dgsfit1, linear_basis

    from srxraylib.plot.gol import plot, plot_table

    # loads file with data to fit
    input_array = numpy.loadtxt("aps_axo_influence_functions2019.dat")

    abscissas = input_array[:, 0].copy()
    print("abscisas: ", abscissas)

    tmp = numpy.loadtxt(filein)
    print(">>>>>>>>>>>>>>>>>>", tmp.shape)
    plot(tmp[:, 0], tmp[:, 1], title="data to fit")
    u = numpy.interp(abscissas, 1000 * tmp[:, 0], tmp[:, 1])
    plot(abscissas, u, title="Result of fit")

    sigma = (abscissas[-1] - abscissas[0])
    g = 15 * numpy.exp(- abscissas ** 2 / 2 / sigma)
    mask = None  # g


    if calculate:
        # prepare input format for orthonormalize_a
        col19 = input_array[:, 0].copy() * 0 + 1
        col20 = numpy.linspace(-1,1,input_array.shape[0])

        a = []
        a.append({'a': col19, 'total_squared': 0})
        a.append({'a': col20, 'total_squared': 0})
        for i in [9, 10, 8, 11, 7, 12, 6, 13, 5, 14, 4, 15, 3, 16, 2, 17, 1, 18]:
            a.append({'a': input_array[:, i], 'total_squared':0})

        plot_table(abscissas, input_array[:, 1:].T, title="influence functions")






        # compute the basis
        b, matrix = orthonormalize_a(a, mask=mask)

        # plot basis
        b_array = numpy.zeros((input_array.shape[0],20))


        for i in range(20):
            b_array[:,i] = b[i]["a"]
        plot_table(abscissas, b_array.T, title="basis functions")


        numpy.savetxt("aps_axo_orthonormal_functions2019.dat",b_array)
        print("File written to disk aps_axo_orthonormal_functions2019.dat")
    else:
        b_array = numpy.loadtxt("aps_axo_orthonormal_functions2019.dat")
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>",b_array.shape)
        b = []
        for i in range(b_array.shape[1]):
            b.append({'a': b_array[:, i], 'total_squared':(b_array[:, i]**2).sum()})

    # perform the fit
    v = linear_2dgsfit1(u, b, mask=mask)
    print("coefficients: ",v)

    # evaluate the fitted data form coefficients and basis
    y = linear_basis(v, b)


    plot(abscissas,u,abscissas,y,legend=["Data","Fit"])

    if fileout != "":
        f = open(fileout,'w')
        for i in range(abscissas.size):
            f.write("%g  %g \n"%(1e-3*abscissas[i],y[i]))
        f.close()
        print("File %s written to disk"%fileout)

    return v


if __name__ == "__main__":
    # filein = "/home/manuel/Oasys/correction.dat"
    # fileout = "/home/manuel/Oasys/correction_fitted.dat"

    root="/home/manuel/OASYS1.2/alsu-scripts/paper-wofry/correctioncryogenic"

    ROOT = ["/home/manuel/OASYS1.2/alsu-scripts/paper-wofry/correctioncryogenic",
            "/home/manuel/OASYS1.2/alsu-scripts/paper-wofry/correctioncryogenicKh",
            "/home/manuel/OASYS1.2/alsu-scripts/paper-wofry/correctionwater1",
            "/home/manuel/OASYS1.2/alsu-scripts/paper-wofry/correctionwater2"]

    for root in ROOT:
        filein = root+".txt"
        fileout = root+"fit.txt"

        v = fit_correction(filein,fileout=fileout,calculate=False)

        print("Coefficients of the orthonormal basis: ")
        v_labels = []
        for i in range(v.size):
            v_labels.append("v[%d]"%i)
            print("v[%d] = %5.2f nm"%(i,1e9*v[i]))

        # import matplotlib.pyplot as plt
        # fig = plt.figure()
        # ax = fig.add_axes([0,0,1,1])
        # ax.bar(v_labels,v)
        # plt.show()
        #
        import matplotlib.pyplot as plt; plt.rcdefaults()
        import numpy as np
        import matplotlib.pyplot as plt

        # objects =   v_labels
        y_pos = np.arange(v.size)


        plt.bar(y_pos, v, align='center', alpha=0.5)
        plt.xticks(y_pos, v_labels)
        plt.ylabel('Usage')
        plt.title('Corfficients in ')

        plt.show()
