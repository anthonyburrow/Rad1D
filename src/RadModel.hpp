#pragma once

#include <string>
#include <vector>
#include <pybind11/pybind11.h>

#include "io.hpp"

namespace myLib
{
    class RadModel
    {
    public:
        radParams params;
        std::vector<feature> lineList;

        std::vector<double> quadMu;
        std::vector<double> quadW;

        std::vector<double> tauCont;
        std::vector<double> T;

        // Methods
        RadModel(const pybind11::dict &dictParams);

        std::vector<std::vector<double>> genSpectrum(bool normalize = false);
        std::vector<std::vector<double>> convergenceTest(const double &lam);

        // Properties
        std::vector<double> getTau();
        std::vector<double> getT();

    private:

    };
}