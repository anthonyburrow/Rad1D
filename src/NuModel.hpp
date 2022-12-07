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

        // Lambda iteration vectors
        std::vector<std::vector<double>> lambda;
        std::vector<double> lambdaA;
        std::vector<double> lambdaB;
        std::vector<double> lambdaC;

        // Methods
        NuModel(const double &lam, const RadModel &radModel);

        double calcFlux();
        double calcF0();
        void iterate(const bool accelerate = true);
        std::vector<std::vector<double>> converge(const bool checkConverged = true,
                                                  const bool returnResults = false);
        std::vector<std::vector<double>> NgConverge(const bool checkConverged = true,
                                                    const bool returnResults = false);

    private:
        void calcTau();
        void setInitialCond();
        void initLambda();
        std::vector<double> SoverB();
    };
}