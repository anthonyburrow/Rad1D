#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "RadModel.hpp"
#include "constants.hpp"
#include "blackbody.hpp"
#include "io.hpp"
#include "util.hpp"
#include "gaussianQuadrature.hpp"
#include "lineProfiles.hpp"

using namespace std;

namespace myLib
{
    vector<vector<double>> initSpectrum(const RadModel &radModel)
    {
        const double &nZones = radModel.params.nZones;
        const double &waveStart = radModel.params.waveStart;
        const double &waveEnd = radModel.params.waveEnd;

        // Add continuum points
        int nWaveCont = int(radModel.params.contRes * (waveEnd - waveStart));
        vector<double> wave(nWaveCont);

        const double slope = (waveEnd - waveStart) / (nWaveCont - 1);
        for (int i = 0; i < nWaveCont; i++)
        {
            wave[i] = i * slope + waveStart;
        }

        // Add line points
        double sigma, lineStart, lineEnd;
        int nWaveLine;
        for (feature line : radModel.lineList)
        {
            // Choose highest temperature to ensure it's resolved
            sigma = gaussianWidth(radModel.T[nZones - 1], line);

            cout << sigma << endl;

            // Go to 3-sigma on each side
            lineStart = line.resonanceWave - 3.0 * sigma;
            lineEnd = line.resonanceWave + 3.0 * sigma;
            nWaveLine = int(radModel.params.lineRes * (waveEnd - waveStart));

            const double slope = (lineEnd - lineStart) / (nWaveLine - 1);
            for (int i = 0; i < nWaveLine; i++)
            {
                wave.push_back(i * slope + lineStart);
            }
        }

        // Setup 2D spectrum vector
        sort(wave.begin(), wave.end());

        const int nWave = wave.size();
        vector<vector<double>> spectrum(nWave, vector<double>(2));

        for (int i = 0; i < nWave; i++)
        {
            spectrum[i][0] = wave[i];
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
