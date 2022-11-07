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
        std::vector<double> tau;
        std::vector<double> B;
        std::vector<double> S;
        std::vector<double> J;
        std::vector<std::vector<double>> lambda;

        // Methods
        NuModel(const double &lam, const RadModel &radModel);

        double calcFlux();
        double calcF0();
        void iterate();

    private:
        void calcTau();
        void setInitialCond();
    };
}