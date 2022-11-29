#include <cmath>
#include <iostream>

#include "constants.hpp"
#include "util.hpp"

using namespace std;

namespace myLib
{
    const double calcHopfQ(const double &tau)
    {
        const double hopfA = 0.7104520194746458;
        const double hopfB = -0.24241757247081414;
        const double hopfC = 0.2215866721657009;

        const vector<double> En = expn(3, tau);
        const double q = hopfA + hopfB * En[2] + hopfC * En[3];

        return q;
    }

    const double calcTemperature(const double &tau, const double &Teff)
    {
        // Look into q = q(tau)
        // const double q = 1.0 / sqrt(3.0);

        const double q = calcHopfQ(tau);
        // cout << tau << " " << q << " " << tau + q << endl;
        const double T = Teff * pow(0.75 * (tau + q), 0.25);

        return T;
    }

    const double planck(const double &lam, const double &T)
    {
        const double bbScale = 2.0 * hc * c * 1e8;
        const double bb = bbScale / (pow(lam, 5) * (exp(1. / (k_hc * T * lam)) - 1));
        return bb;
    }
}