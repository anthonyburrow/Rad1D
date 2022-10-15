#pragma once

#include <vector>

#include "RadModel.hpp"

namespace myLib
{
    void initWave(RadModel &radModel);
    void initT(RadModel &radModel);
    void initTau(RadModel &radModel);

    void initialize(RadModel &radModel);
}