#pragma once

#include <vector>

#include "NuModel.hpp"

namespace myLib
{
    void setupFormalSoln(std::vector<double> &Dtau,
                         std::vector<double> &expDtau,
                         std::vector<double> &e0,
                         std::vector<double> &e1, 
                         std::vector<double> &e2,
                         const double &mu,
                         const std::vector<double> &tau,
                         const int &nZones);

    void calcIpMatrix(std::vector<std::vector<double>> &Ip,
                      const std::vector<double> &Dtau,
                      const std::vector<double> &expDtau,
                      const std::vector<double> &e0,
                      const std::vector<double> &e1, 
                      const std::vector<double> &e2,
                      const int &nZones);

    void calcImMatrix(std::vector<std::vector<double>> &Im,
                      const std::vector<double> &Dtau,
                      const std::vector<double> &expDtau,
                      const std::vector<double> &e0,
                      const std::vector<double> &e1, 
                      const std::vector<double> &e2,
                      const int &nZones);

    void calcLambda(NuModel &nuModel);
}
