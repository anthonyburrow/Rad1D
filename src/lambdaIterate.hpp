#pragma once

#include <vector>

#include "NuModel.hpp"

namespace myLib
{
    void calcJ(NuModel &nuModel);
    void ALIcalcJ(NuModel &nuModel);
    void calcS(NuModel &nuModel);

    void lambdaIteration(NuModel &nuModel);
    void ALI(NuModel &nuModel);
    void NgIteration(NuModel &nuModel,
                     const std::vector<double> &S3,
                     const std::vector<double> &S2,
                     const std::vector<double> &S1,
                     const std::vector<double> &S0);
}