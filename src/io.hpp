#pragma once

#include <string>
#include <vector>
#include <pybind11/pybind11.h>

namespace myLib
{
    // Prevent compiler error when using RadModel types here
    // class RadModel;

    struct radParams
    {
        std::string dataDir = "./data";
        std::string outFilename = "./synthetic.dat";
        double waveStart = 4000.0;
        double waveEnd = 7000.0;
        double tauMax = 100.0;
        double eps = 0.001;
        double Teff = 6000.0;
        int nWave = 1000;
        int nZones = 256;
        int maxIter = 100;
        int nQuad = 8;
    };

    radParams readParams(const std::string &fileName);
    radParams dictToParams(const pybind11::dict &dictParams);
    void printParams(const radParams &params);
}
