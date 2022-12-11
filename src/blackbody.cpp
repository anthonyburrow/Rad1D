#include <cmath>
#include <iostream>

#include "constants.hpp"
#include "util.hpp"

using namespace std;

namespace myLib
{
    const double calcHopfQ(const double &tau)
    {
        const vector<double> En = expn(3, tau);
        const double q = hopfA + hopfB * En[2] + hopfC * En[3];

        return q;
    }

    const double calcTemperature(const double &tau, const double &Teff)
    {
        const double q = calcHopfQ(tau);
        const double T = Teff * pow(0.75 * (tau + q), 0.25);

        return T;
    }

    const double planck(const double &lam, const double &T)
    {
        const double bbScale = 2.0 * hc * c * 1e32;
        const double bb = bbScale / (pow(lam, 5) * (exp(1. / (k_hc * T * lam)) - 1));
        return bb;
    }
}