#pragma once

#include <vector>

#include "RadModel.hpp"

namespace myLib
{
    void calcLambda(const RadModel &radModel, std::vector<double> &A,
                    std::vector<double> &B, std::vector<double> &C);

    void calcLambda(const RadModel &radModel,
                    std::vector<std::vector<double>> &lambda);
}
