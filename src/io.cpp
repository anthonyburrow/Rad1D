#include <iostream>
#include <string>
#include <pybind11/pybind11.h>

#include "io.hpp"

using namespace std;

namespace myLib
{
    radParams dictToParams(const pybind11::dict &dictParams)
    {
        radParams params;

        if (dictParams.contains("data_dir")) {
            params.dataDir = dictParams["data_dir"].cast<string>();
        }
        if (dictParams.contains("wave_start")) {
            params.waveStart = dictParams["wave_start"].cast<double>();
        }
        if (dictParams.contains("wave_end")) {
            params.waveEnd = dictParams["wave_end"].cast<double>();
        }
        if (dictParams.contains("tau_max")) {
            params.tauMax = dictParams["tau_max"].cast<double>();
        }
        if (dictParams.contains("eps")) {
            params.eps = dictParams["eps"].cast<double>();
        }
        if (dictParams.contains("T_eff")) {
            params.Teff = dictParams["T_eff"].cast<double>();
        }
        if (dictParams.contains("n_wave")) {
            params.nWave = dictParams["n_wave"].cast<int>();
        }
        if (dictParams.contains("n_zones")) {
            params.nZones = dictParams["n_zones"].cast<int>();
        }
        if (dictParams.contains("max_iter")) {
            params.maxIter = dictParams["max_iter"].cast<int>();
        }
        if (dictParams.contains("n_quad")) {
            params.nQuad = dictParams["n_quad"].cast<int>();
        }

        printParams(params);

        return params;
    }

    void printParams(const radParams &params)
    {
        cout << "  Data directory:          " << params.dataDir
    <<  endl << "  Wavelength range:        " << params.waveStart << " - " << params.waveEnd
    <<  endl << "  Wavelength points:       " << params.nWave
    <<  endl << "  Thermalization:          " << params.eps
    <<  endl << "  T_eff:                   " << params.Teff
    <<  endl << "  Max tau:                 " << params.tauMax
    <<  endl << "  Number of tau zones:     " << params.nZones
    <<  endl << "  Maximum iterations:      " << params.maxIter
    <<  endl << "  Order of Gaussian quad.: " << params.nQuad
    <<  endl;
    }
}