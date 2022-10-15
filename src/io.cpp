#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
// #include <limits>
// #include <iomanip>
#include <pybind11/pybind11.h>

#include "io.hpp"
#include "constants.hpp"

using namespace std;

namespace myLib
{
    radParams readParams(const string &fileName)
    {
        cout << "Reading from parameter file: " << fileName << endl;

        ifstream paramFile(fileName);
        string line;
        radParams params;

        int count = 0;
        while (getline(paramFile, line))
        {
            if (line[0] == '#' || line[0] == '\0') { continue; }

            stringstream iss(line);

            switch(count)
            {
                case 0 :
                {
                    params.dataDir = line.substr(31);
                    break;
                }
                case 1 :
                {
                    std::string outFilename;
                    iss >> outFilename;
                    // params.outFilename = outFilename;
                    params.outFilename = line;
                    break;
                }
                case 2 :
                    double waveStart, waveEnd;
                    iss >> waveStart >> waveEnd;
                    params.waveStart = waveStart;
                    params.waveEnd = waveEnd;
                    break;
                case 3 :
                    int nWave;
                    iss >> nWave;
                    params.nWave = nWave;
                    break;
                case 4 :
                    double tauMax;
                    iss >> tauMax;
                    params.tauMax = tauMax;
                    break;
                case 5 :
                    int nZones;
                    iss >> nZones;
                    params.nZones = nZones;
                    break;
                case 6 :
                    int maxIter;
                    iss >> maxIter;
                    params.maxIter = maxIter;
                    break;
                case 7 :
                    int nQuad;
                    iss >> nQuad;
                    params.nQuad = nQuad;
                    break;
                case 8 :
                    double eps;
                    iss >> eps;
                    params.eps = eps;
                    break;
                case 9 :
                    double Teff;
                    iss >> Teff;
                    params.Teff = Teff;
                    break;
                default :
                    cout << "Too many lines in param file" << endl;
            }

            count++;
        }
        paramFile.close();

        printParams(params);

        return params;
    }

    radParams dictToParams(const pybind11::dict &dictParams)
    {
        radParams params;

        if (dictParams.contains("data_dir")) {
            params.dataDir = dictParams["data_dir"].cast<string>();
        }
        if (dictParams.contains("out_filename")) {
            params.outFilename = dictParams["out_filename"].cast<string>();
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
            params.eps = dictParams["T_eff"].cast<double>();
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
    <<  endl << "  Output spectrum file:    " << params.outFilename
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