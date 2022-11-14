#pragma once

#include <string>
#include <pybind11/pybind11.h>

// #include "RadModel.hpp"

namespace myLib
{
    // Prevent compiler error when using RadModel types here
    class RadModel;

    struct feature
    {
        double tauRatio;
        double resonanceWave;
        double mass;
    };

    struct radParams
    {
        std::string dataDir = "../data";
        double waveStart = 4000.0;
        double waveEnd = 7000.0;
        double tauMax = 1e6;
        double eps = 1e-4;
        double Teff = 6000.0;
        double contRes = 0.05;
        double lineRes = 3.0;
        double epsConverge = 1e-6;
        int nZones = 256;
        int maxIter = 100;
        int nQuad = 8;
        bool verbose = true;
    };

    radParams dictToParams(const pybind11::dict &dictParams);
    void printParams(const radParams &params);
    void readLineList(RadModel &radModel);
}
