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

        // Physical quantites
        std::vector<double> B;
        std::vector<double> S;
        std::vector<double> J;

        // Methods
        NuModel(const double &lam, const RadModel &radModel);

        double calcFlux();
        double calcF0();
        void iterate();

    private:
        void setBoundary();

    };
}