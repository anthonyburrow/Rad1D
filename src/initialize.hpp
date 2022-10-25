#pragma once

#include <vector>

#include "RadModel.hpp"

namespace myLib
{
    std::vector<std::vector<double>> initSpectrum(const RadModel &radModel);

    void initT(RadModel &radModel);
    void initTau(RadModel &radModel);
    void initialize(RadModel &radModel);
}