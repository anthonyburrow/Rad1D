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
                     std::vector<double> &S3,
                     std::vector<double> &S2,
                     std::vector<double> &S1,
                     std::vector<double> &S0);
}