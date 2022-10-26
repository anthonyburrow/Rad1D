#include <cmath>
#include <iostream>

#include "constants.hpp"

using namespace std;

namespace myLib
{
    const double calcTemperature(const double &tau, const double &Teff)
    {
        // Look into q = q(tau)
        const double q = 1.0 / sqrt(3.0);
        const double T = Teff * pow(0.75 * (tau + q), 0.25);

        return T;
    }

    const double planck(const double &lam, const double &T)
    {
        const double bb = 1. / (pow(lam, 5) * (exp(1. / (k_hc * T * lam)) - 1));
        return bb;
    }
}