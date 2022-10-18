#pragma once

#include <vector>

#include "RadModel.hpp"

namespace myLib
{
    void calcLambda(const RadModel &radModel,
                    std::vector<std::vector<double>> &lambda);
}
