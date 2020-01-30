import numpy

#
# script to make the calculations (created by XOPPY:undulator_radiation)
#
from orangecontrib.xoppy.util.xoppy_undulators import xoppy_calc_undulator_radiation

def run_undulator_radiation(electron_energy=2.0):
    h5_parameters = dict()
    h5_parameters["ELECTRONENERGY"]          = electron_energy
    h5_parameters["ELECTRONENERGYSPREAD"]    = 0.00095
    h5_parameters["ELECTRONCURRENT"]         = 0.5
    h5_parameters["ELECTRONBEAMSIZEH"]       = 0.0002491496538227577
    h5_parameters["ELECTRONBEAMSIZEV"]       = 8.165414870047938e-06
    h5_parameters["ELECTRONBEAMDIVERGENCEH"] = 9.595088728069715e-06
    h5_parameters["ELECTRONBEAMDIVERGENCEV"] = 4.803185217675258e-06
    h5_parameters["PERIODID"]                = 0.015
    h5_parameters["NPERIODS"]                = 133.0
    h5_parameters["KV"]                      = 1.479864
    h5_parameters["KH"]                      = 0.0
    h5_parameters["KPHASE"]                  = 0.0
    h5_parameters["DISTANCE"]                = 20.63
    h5_parameters["SETRESONANCE"]            = 0
    h5_parameters["HARMONICNUMBER"]          = 1
    h5_parameters["GAPH"]                    = 0.0018
    h5_parameters["GAPV"]                    = 0.0018
    h5_parameters["HSLITPOINTS"]             = 41
    h5_parameters["VSLITPOINTS"]             = 41
    h5_parameters["METHOD"]                  = 2
    h5_parameters["PHOTONENERGYMIN"]         = 12000.0
    h5_parameters["PHOTONENERGYMAX"]         = 12001.0
    h5_parameters["PHOTONENERGYPOINTS"]      = 2
    h5_parameters["USEEMITTANCES"]           = 1

    e, h, v, p, code = xoppy_calc_undulator_radiation(
            ELECTRONENERGY           = h5_parameters["ELECTRONENERGY"]         ,
            ELECTRONENERGYSPREAD     = h5_parameters["ELECTRONENERGYSPREAD"]   ,
            ELECTRONCURRENT          = h5_parameters["ELECTRONCURRENT"]        ,
            ELECTRONBEAMSIZEH        = h5_parameters["ELECTRONBEAMSIZEH"]      ,
            ELECTRONBEAMSIZEV        = h5_parameters["ELECTRONBEAMSIZEV"]      ,
            ELECTRONBEAMDIVERGENCEH  = h5_parameters["ELECTRONBEAMDIVERGENCEH"],
            ELECTRONBEAMDIVERGENCEV  = h5_parameters["ELECTRONBEAMDIVERGENCEV"],
            PERIODID                 = h5_parameters["PERIODID"]               ,
            NPERIODS                 = h5_parameters["NPERIODS"]               ,
            KV                       = h5_parameters["KV"]                     ,
            KH                       = h5_parameters["KH"]                     ,
            KPHASE                   = h5_parameters["KPHASE"]                 ,
            DISTANCE                 = h5_parameters["DISTANCE"]               ,
            SETRESONANCE             = h5_parameters["SETRESONANCE"]           ,
            HARMONICNUMBER           = h5_parameters["HARMONICNUMBER"]         ,
            GAPH                     = h5_parameters["GAPH"]                   ,
            GAPV                     = h5_parameters["GAPV"]                   ,
            HSLITPOINTS              = h5_parameters["HSLITPOINTS"]            ,
            VSLITPOINTS              = h5_parameters["VSLITPOINTS"]            ,
            METHOD                   = h5_parameters["METHOD"]                 ,
            PHOTONENERGYMIN          = h5_parameters["PHOTONENERGYMIN"]        ,
            PHOTONENERGYMAX          = h5_parameters["PHOTONENERGYMAX"]        ,
            PHOTONENERGYPOINTS       = h5_parameters["PHOTONENERGYPOINTS"]     ,
            USEEMITTANCES            = h5_parameters["USEEMITTANCES"]          ,
            h5_file                  = "undulator_radiation.h5",
            h5_entry_name            = "XOPPY_RADIATION",
            h5_initialize            = True,
            h5_parameters            = h5_parameters,
            )
    return p,h,v

if __name__ == "__main__":
    # example plot
    from srxraylib.plot.gol import plot_image


    electron_energy = 1.902
    p,h,v = run_undulator_radiation(electron_energy=electron_energy)
    plot_image(p[0],h,v,title="Flux [photons/s] per 0.1 bw per mm2 at %9.3f eV"%(12000.0),xtitle="H [mm]",ytitle="V [mm]")
    #
    # end script
#
