#pragma once

#include <string>
#include <vector>

namespace myLib
{
    const double trapezoidal(const std::vector<double> X,
                             const std::vector<double> integrand);

    const int closestIndex(const std::vector<double> &vec,
                           const double value);
}