#include <iostream>
#include <cmath>
#include <vector>

#include "lambda.hpp"
#include "NuModel.hpp"
#include "util.hpp"
#include "constants.hpp"

using namespace std;

namespace myLib
{
    void setupFormalSoln(vector<double> &Dtau, vector<double> &expDtau,
                         vector<double> &e0, vector<double> &e1, 
                         vector<double> &e2,
                         const double &mu, const vector<double> &tau,
                         const int &nZones)
    {
        for (int i = 1; i < nZones; i++)
        {
            Dtau[i - 1] = (tau[i] - tau[i - 1]) / mu;

            if (Dtau[i - 1] < 1e-8)
            {
                expDtau[i - 1] = 1.0 - Dtau[i - 1];

                e0[i] = Dtau[i - 1];
                e1[i] = zero;
                e2[i] = pow(Dtau[i - 1], 2.0);
            }
            else
            {
                expDtau[i - 1] = exp(-Dtau[i - 1]);

                e0[i] = 1.0 - expDtau[i - 1];
                e1[i] = Dtau[i - 1] - e0[i];
                e2[i] = pow(Dtau[i - 1], 2.0) - 2.0 * e1[i];
            }
        }
    }

    void calcIpMatrix(vector<vector<double>> &Ip,
                      const vector<double> &Dtau, const vector<double> &expDtau,
                      const vector<double> &e0, const vector<double> &e1, 
                      const vector<double> &e2, const int &nZones)
    {
        vector<double> alphaP(nZones, zero);
        vector<double> betaP(nZones, zero);
        vector<double> gammaP(nZones, zero);

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

        // Use linear interpolation for columns 1 and N
        alphaP[nZones - 1] = zero;
        betaP[0] = e1[1] / Dtau[0];
        betaP[nZones - 1] = 1.0;   // I^+(tau_max) = B(T(tau_max))
        gammaP[0] = e0[1] - e1[1] / Dtau[0];

        // Correct for under-correction in interpolation:
        // (These correspond to I+- values so they cannot be < 0)
        for (int i = 0; i < nZones; i++)
        {
            // alphaP[i] = max(alphaP[i], zero);
            // betaP[i] = max(betaP[i], zero);
            // gammaP[i] = max(gammaP[i], zero);
        }

        for (int i = 1; i < nZones - 1; i++)
        {
            // A (i + 1, i)
            Ip[i + 1][i] = alphaP[i + 1];

            // B (i, i)
            Ip[i][i] = Ip[i + 1][i] * expDtau[i] + betaP[i];

            // C (i - 1, i)
            Ip[i - 1][i] = Ip[i][i] * expDtau[i - 1] + gammaP[i - 1];

            // I^+ from row (i - 2) to 1
            for (int k = i - 2; k > -1; k--)
            {
                Ip[k][i] = Ip[k + 1][i] * expDtau[k];
            }
        }

        // Column 1 - Use linear interpolation

        Ip[0][0] = betaP[0];

        // Column N - Use linear interpolation

        Ip[nZones - 1][nZones - 1] = betaP[nZones - 1];

        Ip[nZones - 2][nZones - 1] = Ip[nZones - 1][nZones - 1] * expDtau[nZones - 2] +
                                     gammaP[nZones - 2];

        for (int k = nZones - 3; k > -1; k--)
        {
            Ip[k][nZones - 1] = Ip[k + 1][nZones - 1] * expDtau[k];
        }
    }

    void calcImMatrix(vector<vector<double>> &Im,
                      const vector<double> &Dtau, const vector<double> &expDtau,
                      const vector<double> &e0, const vector<double> &e1, 
                      const vector<double> &e2, const int &nZones)
    {
        vector<double> alphaM(nZones, zero);
        vector<double> betaM(nZones, zero);
        vector<double> gammaM(nZones, zero);

        for (int i = 1; i < nZones - 1; i++)
        {
            alphaM[i] = e0[i] +
                        (e2[i] - (Dtau[i] + 2 * Dtau[i - 1]) * e1[i]) /
                        (Dtau[i - 1] * (Dtau[i] + Dtau[i - 1]));
            betaM[i] = ((Dtau[i] + Dtau[i - 1]) * e1[i] - e2[i]) /
                        (Dtau[i] * Dtau[i - 1]);
            gammaM[i] = (e2[i] - Dtau[i - 1] * e1[i]) /
                        (Dtau[i] * (Dtau[i] + Dtau[i - 1]));
        }

        // Use linear interpolation for columns 1 and N
        alphaM[nZones - 1] = e0[nZones - 1] - e1[nZones - 1] / Dtau[nZones - 2];   
        betaM[nZones - 1] = e1[nZones - 1] / Dtau[nZones - 2];
        gammaM[0] = zero;

        // Correct for under-correction in interpolation:
        // (These correspond to I+- values so they cannot be < 0)
        for (int i = 0; i < nZones; i++)
        {
            // alphaM[i] = max(alphaM[i], zero);
            // betaM[i] = max(betaM[i], zero);
            // gammaM[i] = max(gammaM[i], zero);
        }

        // Columns 2 through (N - 1)
        for (int i = 1; i < nZones - 1; i++)
        {
            // C (i - 1, i)
            Im[i - 1][i] = gammaM[i - 1];

            // B (i, i)
            Im[i][i] = Im[i - 1][i] * expDtau[i - 1] + betaM[i];

            // A (i + 1, i)
            Im[i + 1][i] = Im[i][i] * expDtau[i] + alphaM[i + 1];

            // I^- from row (i + 2) to N
            for (int k = i + 2; k < nZones; k++)
            {
                Im[k][i] = Im[k - 1][i] * expDtau[k - 1];
            }
        }

        // Column 1 - Use linear interpolation

        // i^-(0) = 0

        // Use linear alphaM[1] instead of the parabolic one?
        Im[1][0] = e0[1] - e1[1] / Dtau[0];
        // Im[1][0] = alphaM[1];

        for (int k = 2; k < nZones; k++)
        {
            Im[k][0] = Im[k - 1][0] * expDtau[k - 1];
        }

        // Column N - Use linear interpolation

        Im[nZones - 1][nZones - 1] = betaM[nZones - 1];
    }

    void calcLambda(const NuModel &nuModel, vector<vector<double>> &lambda)
    {
        const int &nZones = nuModel.params.nZones;
        const int &nQuad = nuModel.params.nQuad;
        const int halfQuad = int(0.5 * nQuad);
        const vector<double> &tau = nuModel.tau;
        const vector<double> &quadMu = nuModel.radModel.quadMu;
        const vector<double> &quadW = nuModel.radModel.quadW;

        vector<double> Dtau(nZones, zero);
        vector<double> expDtau(nZones, zero);

        vector<double> e0(nZones, zero);
        vector<double> e1(nZones, zero);
        vector<double> e2(nZones, zero);

        vector<vector<double>> Ip(nZones, vector<double>(nZones, zero));
        vector<vector<double>> Im(nZones, vector<double>(nZones, zero));

        for (int j = 0; j < halfQuad; j++)
        {
            // Setup values needed for both I+ and I-
            setupFormalSoln(Dtau, expDtau, e0, e1, e2,
                            quadMu[j], tau, nZones);

            // Solve for I+ and I- parts of Lambda
            calcIpMatrix(Ip, Dtau, expDtau, e0, e1, e2, nZones);
            calcImMatrix(Im, Dtau, expDtau, e0, e1, e2, nZones);

            // Integrate mu (Quadrature sum) to get J
            for (int i = 0; i < nZones; i++)
            {
                for (int k = 0; k < nZones; k++)
                {
                    lambda[i][k] += quadW[j] * (Im[i][k] + Ip[i][k]);
                }
            }
        }

        for (int i=0; i < nZones; i++)
        {
            for (int j=0; j < nZones; j++)
            {
                lambda[i][j] *= 0.5;
            }
        }
    }
}