#include <iostream>
#include <vector>
#include <cmath>

#include "RadModel.hpp"
#include "constants.hpp"
#include "blackbody.hpp"
#include "io.hpp"
#include "gaussianQuadrature.hpp"

using namespace std;

namespace myLib
{
    void initWave(RadModel &radModel)
    {
        // Uniform wavelengths (could select around line list later)
        const double &waveStart = radModel.params.waveStart;
        const double &waveEnd = radModel.params.waveEnd;
        const int &nWave = radModel.params.nWave;

        const double slope = (waveEnd - waveStart) / (nWave - 1);

        for (int i = 0; i < nWave; i++)
        {
            radModel.spectrum[i][0] = i * slope + waveStart;
        }
    }

    void initTau(RadModel &radModel)
    {
        const int &nZones = radModel.params.nZones;
        const double logTauMin = -7.0;
        const double logTauMax = log10(radModel.params.tauMax);
        const double expSlope = (logTauMax - logTauMin) / (nZones - 2);

        radModel.tau[0] = zero;

        for (int i = 1; i < nZones; i++)
        {
            radModel.tau[i] = (i - 1) * expSlope + logTauMin;
            radModel.tau[i] = pow(10, radModel.tau[i]);
        }
    }

    void initT(RadModel &radModel)
    {
        const int &nZones = radModel.params.nZones;
        const double &Teff = radModel.params.Teff;

        // Look into q = q(tau)
        const double q = 1.0 / sqrt(3.0);
        for (int i = 0; i < nZones; i++)
        {
            radModel.T[i] = Teff * pow(0.75 * (radModel.tau[i] + q), 0.25);
        }
    }

    void initialize(RadModel &radModel)
    {
        initWave(radModel);
        initTau(radModel);

        initT(radModel);

        getWeights(radModel);
    }
}
