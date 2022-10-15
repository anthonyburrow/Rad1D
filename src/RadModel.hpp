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

        // Physical arrays
        std::vector<std::vector<double>> spectrum;
        std::vector<double> tau;
        std::vector<double> T;

        // Methods
        RadModel(const radParams &params);
        RadModel(const pybind11::dict &dictParams);

        std::vector<std::vector<double>> genSpectrum();
        void rescaleFlux();

        // Properties
        std::vector<double> getTau();
        std::vector<double> getT();

    private:

    };
}