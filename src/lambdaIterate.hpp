#pragma once

#include <vector>

#include "NuModel.hpp"

namespace myLib
{
    void calcJ(NuModel &nuModel);
    void calcS(NuModel &nuModel);

    void lambdaIteration(NuModel &nuModel);
}