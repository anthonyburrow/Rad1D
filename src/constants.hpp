#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace myLib
{
    // Math constants
    // static const double zero = 1e-99;
    static const double zero = 0.0;
    static const double &pi = M_PI;
    static const double pi4 = 4.0 * M_PI;
    static const double pi4_sq = pi4 * pi4;
    static const double one_third = 1.0 / 3.0;
    static const double four_thirds = 4.0 / 3.0;
    static const double pi4_3 = pi4 * one_third;
    static const double inv_sqrt_pi = 1.0 / sqrt(pi);

    // Physical constants
    static const double hc = 1.98644586e-16;               // erg cm
    static const double hc_eV = 12398.4198;                // eV A
    static const double k = 8.617333262e-5;                // eV K^{-1}
    static const double k_hc = k / hc_eV;                  // K^{-1} A^{-1}
    static const double c = 2.99792458e10;                 // cm s^{-1}
    static const double constDB = 4.301415e-7;             // K^{-1/2} amu^{1/2}

    static const double hopfA = 0.7104520194746458;
    static const double hopfB = -0.24241757247081414;
    static const double hopfC = 0.2215866721657009;

    // Settings
    static const double expDtauThreshold = 1e-7;      // Dtau where exponentials are expanded
    static const double interpDtauThreshold = 1e-2;   // Dtau where linear S interpolation is used
}