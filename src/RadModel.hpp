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

        std::vector<double> quadMu;
        std::vector<double> quadW;

        std::vector<std::vector<double>> spectrum;
        std::vector<double> tau;
        std::vector<double> T;

        std::vector<std::vector<double>> lambda;

        // Methods
        RadModel(const pybind11::dict &dictParams);

        void initLambda();
        std::vector<std::vector<double>> genSpectrum(bool normalize = false);
        std::vector<std::vector<double>> convergenceTest(const double &lam);

        // Properties
        std::vector<double> getTau();
        std::vector<double> getT();

    private:
        void normalizeFlux();
    };
}