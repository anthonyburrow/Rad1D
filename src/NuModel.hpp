#pragma once

#include <vector>

#include "RadModel.hpp"
#include "io.hpp"

namespace myLib
{
    class NuModel
    {
    public:
        const double &lam;
        const RadModel &radModel;
        const radParams &params;

        // Physical arrays
        std::vector<double> B;
        std::vector<double> S;
        // std::vector<std::vector<double>> I;
        std::vector<double> J;

        // Methods
        NuModel(const double &lam, const RadModel &radModel);

        const double getFlux();

    private:
        void setBoundary();

    };
}