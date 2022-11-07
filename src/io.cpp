#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <pybind11/pybind11.h>

#include "RadModel.hpp"
#include "io.hpp"

using namespace std;

namespace myLib
{
    radParams dictToParams(const pybind11::dict &dictParams)
    {
        radParams params;

        if (dictParams.contains("verbose")) {
            params.verbose = dictParams["verbose"].cast<bool>();
        }

        if (params.verbose) { cout << "Reading parameters..." << endl; }

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
        if (dictParams.contains("cont_res")) {
            params.contRes = dictParams["cont_res"].cast<double>();
        }
        if (dictParams.contains("line_res")) {
            params.lineRes = dictParams["line_res"].cast<double>();
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
        if (dictParams.contains("eps_converge")) {
            params.epsConverge = dictParams["eps_converge"].cast<double>();
        }

        if (params.verbose) { printParams(params); }

        return params;
    }

    void printParams(const radParams &params)
    {
        cout << "  Data directory:                 " << params.dataDir
    <<  endl << "  Wavelength range [A]:           " << params.waveStart << " - " << params.waveEnd
    <<  endl << "  Continuum resolution [pts / A]: " << params.contRes
    <<  endl << "  Line resolution [pts / A]:      " << params.lineRes
    <<  endl << "  Thermalization:                 " << params.eps
    <<  endl << "  Effective temperature [K]:      " << params.Teff
    <<  endl << "  Max tau:                        " << params.tauMax
    <<  endl << "  Number of tau points:           " << params.nZones
    <<  endl << "  Maximum iterations:             " << params.maxIter
    <<  endl << "  Order of Gaussian quad.:        " << params.nQuad
    <<  endl << "  Epsilon for convergence:        " << params.epsConverge
    <<  endl;
    }

    void readLineList(RadModel &radModel)
    {
        vector<feature> &lineList = radModel.lineList;
        string fileName = radModel.params.dataDir + "/lines.dat";

        ostringstream output;
        output << "Reading from line list: " << fileName;
        radModel.log(output);

        ifstream lineListFile(fileName);
        string line;

        double tauRatio, resonanceWave, mass;

        while (getline(lineListFile, line))
        {
            if (line[0] == '#' || line[0] == '\0') { continue; }

            stringstream iss(line);
            iss >> tauRatio >> resonanceWave >> mass;

            feature newFeature = feature();
            newFeature.tauRatio = tauRatio;
            newFeature.resonanceWave = resonanceWave;
            newFeature.mass = mass;

            lineList.push_back(newFeature);
        }
        lineListFile.close();
    }
}