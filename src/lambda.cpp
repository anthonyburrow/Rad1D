#include <iostream>
#include <cmath>
#include <vector>

#include "lambda.hpp"
#include "RadModel.hpp"
#include "util.hpp"
#include "constants.hpp"

using namespace std;

namespace myLib
{
    void calcLambda(const RadModel &radModel, vector<vector<double>> &lambda)
    {
        const int &nZones = radModel.params.nZones;
        const int &nQuad = radModel.params.nQuad;
        const int halfQuad = int(0.5 * nQuad);
        const vector<double> &tau = radModel.tau;
        const vector<double> &quadMu = radModel.quadMu;
        const vector<double> &quadW = radModel.quadW;
        double Ip;
        double Im;
        double IPrev;

        vector<double> Dtau(nZones, zero);
        vector<double> expDtau(nZones, zero);

        vector<double> e0(nZones, zero);
        vector<double> e1(nZones, zero);
        vector<double> e2(nZones, zero);

        vector<double> alphaP(nZones, zero);
        vector<double> betaP(nZones, zero);
        vector<double> gammaP(nZones, zero);
        vector<double> alphaM(nZones, zero);
        vector<double> betaM(nZones, zero);
        vector<double> gammaM(nZones, zero);


        for (int j = 0; j < halfQuad; j++)
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
                alphaM[i] = e0[i] +
                            (e2[i] - (Dtau[i] + 2 * Dtau[i - 1]) * e1[i]) /
                            (Dtau[i - 1] * (Dtau[i] + Dtau[i - 1]));
                betaM[i] = ((Dtau[i] + Dtau[i - 1]) * e1[i] - e2[i]) /
                           (Dtau[i] * Dtau[i - 1]);
                gammaM[i] = (e2[i] - Dtau[i - 1] * e1[i]) /
                            (Dtau[i] * (Dtau[i] + Dtau[i - 1]));
            }

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

            // tbh I have no idea for these...
            // I just used the linear interpolation values for them
            alphaM[nZones - 1] = e0[nZones - 1] - e1[nZones - 1] / Dtau[nZones - 2];   
            alphaP[nZones - 1] = zero;

            betaP[0] = e1[1] / Dtau[0];
            betaM[nZones - 1] = e1[nZones - 1] / Dtau[nZones - 2];

            betaP[nZones - 1] = 1.0;   // I^+(tau_max) normalized

            gammaM[0] = zero;
            gammaP[0] = e0[1] - e1[1] / Dtau[0];

            // Correct for under-correction in interpolation:
            // (These correspond to I+- values so they cannot be < 0)
            for (int i = 0; i < nZones; i++)
            {
                alphaP[i] = max(alphaP[i], zero);
                betaP[i] = max(betaP[i], zero);
                gammaP[i] = max(gammaP[i], zero);
                alphaM[i] = max(alphaM[i], zero);
                betaM[i] = max(betaM[i], zero);
                gammaM[i] = max(gammaM[i], zero);
            }

            // Integrate mu first (Quadrature) to fill Lambda

            // Columns 2 through (N - 1)
            for (int i = 1; i < nZones - 1; i++)
            {
                // C (i - 1, i)
                Im = gammaM[i - 1];
                Ip = (alphaP[i + 1] * expDtau[i] + betaP[i]) * expDtau[i - 1] +
                     gammaP[i - 1];

                lambda[i - 1][i] += quadW[j] * (Im + Ip);

                // I^+ from row (i - 2) to 1 (I^- = 0)
                IPrev = Ip;
                for (int k = i - 2; k > -1; k--)
                {
                    Ip = IPrev * expDtau[k];
                    lambda[k][i] += quadW[j] * Ip;
                    IPrev = Ip;
                }

                // B (i, i)
                Im = gammaM[i - 1] * expDtau[i - 1] + betaM[i];
                Ip = alphaP[i + 1] * expDtau[i] + betaP[i];

                lambda[i][i] += quadW[j] * (Im + Ip);

                // A (i + 1, i)
                Im = (gammaM[i - 1] * expDtau[i - 1] + betaM[i]) * expDtau[i] +
                     alphaM[i + 1];
                Ip = alphaP[i + 1];

                lambda[i + 1][i] += quadW[j] * (Im + Ip);

                // I^- from row (i + 2) to N (I^+ = 0)
                IPrev = Im;
                for (int k = i + 2; k < nZones; k++)
                {
                    Im = IPrev * expDtau[k - 1];
                    lambda[k][i] += quadW[j] * Im;
                    IPrev = Im;
                }
            }

            // Column 1 - Use linear interpolation

            // i^-(0) = 0 so this is just i^+(0) in integral
            lambda[0][0] += quadW[j] * betaP[0];

            // i^-(0) = 0 so i^-_{1 + 1} = alphaM_{1 + 1}
            // Use linear alphaM[1] and alphaP[1] (= 0)
            Im = e0[1] - e1[1] / Dtau[0];
            // Im = alphaM[1];
            lambda[1][0] += quadW[j] * Im;

            IPrev = Im;
            for (int k = 2; k < nZones; k++)
            {
                Im = IPrev * expDtau[k - 1];
                lambda[k][0] += quadW[j] * Im;
                IPrev = Im;
            }

            // Column N - Use linear interpolation

            Im = betaM[nZones - 1];
            Ip = betaP[nZones - 1];
            lambda[nZones - 1][nZones - 1] += quadW[j] * (Im + Ip);

            Ip = betaP[nZones - 1] * expDtau[nZones - 2] + gammaP[nZones - 2];
            lambda[nZones - 2][nZones - 1] += quadW[j] * Ip;

            IPrev = Ip;
            for (int k = nZones - 3; k > -1; k--)
            {
                Ip = IPrev * expDtau[k];
                lambda[k][nZones - 1] += quadW[j] * Ip;
                IPrev = Ip;
            }
        }

        for (int i=0; i < nZones; i++)
        {
            for (int j=0; j < nZones; j++)
            {
                lambda[i][j] *= 0.5;
            }
        }

        // for (int i=0; i < nZones; i++)
        // {
        //     for (int j=0; j < nZones; j++)
        //     {
        //         cout << lambda[i][j] << "  ";
        //     }
        //     cout << endl;
        // }
    }
}