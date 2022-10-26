#include <iostream>
#include <vector>
#include <cmath>

#include "RadModel.hpp"
#include "constants.hpp"
#include "blackbody.hpp"
#include "io.hpp"
#include "util.hpp"
#include "gaussianQuadrature.hpp"

using namespace std;

namespace myLib
{
    vector<vector<double>> initSpectrum(const RadModel &radModel)
    {
        // Uniform wavelengths (could select around line list later)
        const double &waveStart = radModel.params.waveStart;
        const double &waveEnd = radModel.params.waveEnd;

        // Get number of continuum points (+ line points later)
        const int &nWave = int(radModel.params.contRes * (waveEnd - waveStart));

        vector<vector<double>> spectrum(nWave, vector<double>(2));

        const double slope = (waveEnd - waveStart) / (nWave - 1);

        for (int i = 0; i < nWave; i++)
        {
            spectrum[i][0] = i * slope + waveStart;
        }

        return spectrum;
    }

    void initTau(RadModel &radModel)
    {
        const int &nZones = radModel.params.nZones;
        const double logTauMin = -7.0;
        const double logTauMax = log10(radModel.params.tauMax);
        const double expSlope = (logTauMax - logTauMin) / (nZones - 2);

        radModel.tauCont[0] = zero;

        for (int i = 1; i < nZones; i++)
        {
            radModel.tauCont[i] = (i - 1) * expSlope + logTauMin;
            radModel.tauCont[i] = pow(10, radModel.tauCont[i]);
        }

        // Replace closest tau value to 1 with 1
        const int ind = closestIndex(radModel.tauCont, 1.0);
        radModel.tauCont[ind] = 1.0;
    }

    void initT(RadModel &radModel)
    {
        const int &nZones = radModel.params.nZones;
        const double &Teff = radModel.params.Teff;

        for (int i = 0; i < nZones; i++)
        {
            radModel.T[i] = calcTemperature(radModel.tauCont[i], Teff);
        }
    }

    void initialize(RadModel &radModel)
    {
        initTau(radModel);
        initT(radModel);
        getWeights(radModel);
    }
}
