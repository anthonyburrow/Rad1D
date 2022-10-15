#include <cmath>
#include <iostream>

#include "constants.hpp"

using namespace std;

namespace myLib
{
    const double planck(const double &lam, const double &T)
    {
        const double bb = 1. / (pow(lam, 5) * (exp(1. / (k_hc * T * lam)) - 1));
        return bb;
    }
}