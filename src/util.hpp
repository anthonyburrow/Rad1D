#pragma once

#include <string>
#include <vector>

namespace myLib
{
    const double trapezoidal(const std::vector<double> X,
                             const std::vector<double> integrand);

    void normalizeSpec(std::vector<std::vector<double>> &spectrum);

    const int closestIndex(const std::vector<double> &vec,
                           const double value);

    void TridiagonalSoln(std::vector<double> &y,
                         const std::vector<std::vector<double>> &lambda,
                         const std::vector<double> &x,
                         const double &eps);

    const std::vector<double> expn(const int &n,
                                   const double &x);
}