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
        S(params.nZones, zero),
        J(params.nZones, zero)
    {
        setBoundary();
    }

    void NuModel::setBoundary()
    {
        const int &nZones = params.nZones;
        double bb;

        // Initialize S(tau) = B(T(tau)) = 1
        for (int i = 0; i < nZones; i++)
        {
            S[i] = 1.0;
        }
    }

    void NuModel::iterate()
    {
        lambdaIteration(*this);
    }

    const double NuModel::calcFlux()
    {
        const int &maxIter = params.maxIter;

        // Converge S & J
        for (int i=0; i < maxIter; i++)
        {
            iterate();

            // Check/break for convergence?
        }

        // Calc F based on converged S & J
        return calcF0();
    }

    const double NuModel::calcF0()
    {
        const int &nZones = params.maxIter;
        const int halfQuad = int(0.5 * params.nQuad);
        const vector<double> &tau = radModel.tau;
        const vector<double> &quadMu = radModel.quadMu;
        const vector<double> &quadW = radModel.quadW;
        vector<vector<double>> I(nZones, vector<double>(halfQuad));
        double F0 = zero;

        // Just do this for now, optimize later cuz I'm over it rn
        vector<vector<double>> Dtau(nZones, vector<double>(halfQuad));
        vector<vector<double>> expDtau(nZones, vector<double>(halfQuad));

        vector<vector<double>> e0(nZones, vector<double>(halfQuad));
        vector<vector<double>> e1(nZones, vector<double>(halfQuad));
        vector<vector<double>> e2(nZones, vector<double>(halfQuad));

        vector<vector<double>> alphaP(nZones, vector<double>(halfQuad));
        vector<vector<double>> betaP(nZones, vector<double>(halfQuad));
        vector<vector<double>> gammaP(nZones, vector<double>(halfQuad));

        // Calc delta tau's
        for (int i=0; i < nZones - 1; i++)
        {
            for (int j=0; j < halfQuad; j++)
            {
                Dtau[i][j] = (tau[i + 1] - tau[i]) / quadMu[j];
                expDtau[i][j] = exp(-Dtau[i][j]);
            }
        }

        // Calc e's
        for (int i=1; i < nZones; i++)
        {
            for (int j=0; j < halfQuad; j++)
            {
                e0[i][j] = 1.0 - expDtau[i - 1][j];
                e1[i][j] = Dtau[i - 1][j] - e0[i][j];
                e2[i][j] = pow(Dtau[i - 1][j], 2.0) - 2.0 * e1[i][j];
            }
        }

        // Calc parabolic coefficients (from zone 2 to N - 1)
        for (int i=1; i < nZones - 1; i++)
        {
            for (int j=0; j < halfQuad; j++)
            {
                alphaP[i][j] = (e2[i + 1][j] - Dtau[i][j] * e1[i + 1][j]) /
                               (Dtau[i - 1][j] * (Dtau[i][j] + Dtau[i - 1][j]));
                betaP[i][j] = ((Dtau[i][j] + Dtau[i - 1][j]) * e1[i + 1][j] - e2[i + 1][j]) /
                              ((Dtau[i][j] * Dtau[i - 1][j]));
                gammaP[i][j] = e0[i + 1][j] +
                               (e2[i + 1][j] - (Dtau[i - 1][j] + 2 * Dtau[i][j]) * e1[i + 1][j]) /
                               (Dtau[i][j] * (Dtau[i][j] + Dtau[i - 1][j]));
            }
        }

        // Calculate I^+ from inside out
        for (int j=0; j < halfQuad; j++)
        {
            I[nZones - 1][j] = S[nZones - 1];
        }

        for (int i=nZones - 2; i > 0; i--)
        {
            for (int j=0; j < halfQuad; j++)
            {
                I[i][j] = I[i + 1][j] * expDtau[i][j] +
                          alphaP[i][j] * S[i - 1] + betaP[i][j] * S[i] + gammaP[i][j] * S[i + 1];
            }
        }

        for (int j=0; j < halfQuad; j++)
        {
            I[0][j] = I[1][j] * expDtau[0][j] +
                      (e1[1][j] / Dtau[0][j]) * S[0] + (e0[1][j] - e1[1][j] / Dtau[0][j]) * S[1];
        }

        // F(0) = 4 pi H(0) = 2 pi mu_j I_0j Wj
        for (int j=0; j < halfQuad; j++)
        {
            F0 += quadMu[j] * I[0][j] * quadW[j];
        }

        F0 *= 2.0 * pi;

        return F0;
    }
}
