#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace myLib
{
    // Math constants
    static const double zero = 1e-70;
    static const double &pi = M_PI;
    static const double pi4 = 4.0 * M_PI;
    static const double pi4_sq = pi4 * pi4;
    static const double one_third = 1.0 / 3.0;
    static const double four_thirds = 4.0 / 3.0;
    static const double pi4_3 = pi4 * one_third;

    // Physical constants
    static const double hc = 12398.4198;      // eV A
    static const double k = 8.617333262e-5;   // eV K^{-1}
    static const double k_hc = k / hc;        // K^{-1} A^{-1}
    static const double c = 29979.2458;       // cm s^{-1}
}