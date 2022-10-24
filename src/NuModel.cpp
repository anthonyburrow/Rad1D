#include <iostream>
#include <cmath>

#include "NuModel.hpp"
#include "constants.hpp"
#include "blackbody.hpp"
#include "lambdaIterate.hpp"

using namespace std;

namespace myLib
{
    NuModel::NuModel(const double &lam, const RadModel &radModel)
      : lam(lam),
        radModel(radModel),
        params(radModel.params),
        // Allocate vectors
        B(params.nZones, zero),
        S(params.nZones, zero),
        J(params.nZones, zero)
    {
        setBoundary();
    }

    void NuModel::setBoundary()
    {
        const int &nZones = params.nZones;
        double bb;
        const double bbScale = 2.0 * hc * c * 1e8;

        // Initialize S(tau) = B(T(tau)) = 1
        for (int i = 0; i < nZones; i++)
        {
            B[i] = bbScale * planck(lam, radModel.T[i]);
            S[i] = 1.0;
        }
    }

    void NuModel::iterate()
    {
        lambdaIteration(*this);
    }

    double NuModel::calcFlux()
    {
        const int &maxIter = params.maxIter;
        double flux;

        // Converge S & J
        for (int i=0; i < maxIter; i++)
        {
            iterate();

            // Check/break for convergence?
        }

        // Calc F based on converged S & J
        flux = calcF0();
        flux *= B[0];

        return flux;
    }

    double NuModel::calcF0()
    {
        const int &nZones = params.maxIter;
        const int halfQuad = int(0.5 * params.nQuad);
        const vector<double> &tau = radModel.tau;
        const vector<double> &quadMu = radModel.quadMu;
        const vector<double> &quadW = radModel.quadW;
        double Ip = zero;
        double F0 = zero;

        vector<double> Dtau(nZones, zero);
        vector<double> expDtau(nZones, zero);

        vector<double> e0(nZones, zero);
        vector<double> e1(nZones, zero);
        vector<double> e2(nZones, zero);

        vector<double> alphaP(nZones, zero);
        vector<double> betaP(nZones, zero);
        vector<double> gammaP(nZones, zero);

        for (int j=0; j < halfQuad; j++)
        {
            // Calc delta tau's
            for (int i = 1; i < nZones; i++)
            {
                Dtau[i - 1] = (tau[i] - tau[i - 1]) / quadMu[j];
                expDtau[i - 1] = exp(-Dtau[i - 1]);
            }

            // Calc e's
            for (int i = 1; i < nZones; i++)
            {
                e0[i] = 1.0 - expDtau[i - 1];
                e1[i] = Dtau[i - 1] - e0[i];
                e2[i] = pow(Dtau[i - 1], 2.0) - 2.0 * e1[i];
            }

            // Calc parabolic coefficients (from zone 2 to N - 1)
            for (int i = 1; i < nZones - 1; i++)
            {
                alphaP[i] = (e2[i + 1] - Dtau[i] * e1[i + 1]) /
                            (Dtau[i - 1] * (Dtau[i] + Dtau[i - 1]));
                betaP[i] = ((Dtau[i] + Dtau[i - 1]) * e1[i + 1] - e2[i + 1]) /
                           ((Dtau[i] * Dtau[i - 1]));
                gammaP[i] = e0[i + 1] +
                            (e2[i + 1] - (Dtau[i - 1] + 2 * Dtau[i]) * e1[i + 1]) /
                            (Dtau[i] * (Dtau[i] + Dtau[i - 1]));
            }

            // Use linear interpolation values for bounds
            alphaP[nZones - 1] = zero;
            betaP[0] = e1[1] / Dtau[0];
            betaP[nZones - 1] = 1.0;   // I^+(tau_max) normalized
            gammaP[0] = e0[1] - e1[1] / Dtau[0];

            // Correct for under-correction in interpolation:
            // (These correspond to I+- values so they cannot be < 0)
            for (int i = 0; i < nZones; i++)
            {
                // alphaP[i] = max(alphaP[i], zero);
                // betaP[i] = max(betaP[i], zero);
                // gammaP[i] = max(gammaP[i], zero);
                // alphaM[i] = max(alphaM[i], zero);
                // betaM[i] = max(betaM[i], zero);
                // gammaM[i] = max(gammaM[i], zero);
            }

            // Calculate I^+ from inside out
            Ip = 1.0;

            for (int i=nZones - 2; i > 0; i--)
            {
                Ip = Ip * expDtau[i] +
                     alphaP[i] * S[i - 1] + betaP[i] * S[i] + gammaP[i] * S[i + 1];
            }

            Ip = Ip * expDtau[0] +
                 betaP[0] * S[0] + gammaP[0] * S[1];

            // F(0) = 4 pi H(0) = 2 pi mu_j I_0j Wj
            F0 += quadMu[j] * Ip * quadW[j];
        }

        F0 *= 2.0 * pi;

        // cout << F0 << endl;

        return F0;
    }
}
