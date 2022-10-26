#include <iostream>
#include <cmath>

#include "io.hpp"
#include "constants.hpp"

using namespace std;

namespace myLib
{
    const double gaussianProfile(const double &lam, const double &T,
                                 const feature &line)
    {
        const double &lam0 = line.resonanceWave;
        const double &m = line.mass;

        const double deltaLam = constDB * lam0 * sqrt(T / m);

        const double phi = inv_sqrt_pi * exp(-pow((lam - lam0) / deltaLam, 2.0))
                           / deltaLam;
        
        return phi;
    }
}