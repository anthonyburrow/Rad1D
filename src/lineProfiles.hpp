#pragma once

#include "io.hpp"

namespace myLib
{
    const double gaussianWidth(const double &T,
                               const feature &line);

    const double gaussianProfile(const double &lam,
                                 const double &T,
                                 const feature &line);
}