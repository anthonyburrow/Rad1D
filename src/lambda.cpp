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
    void calcLambda(const RadModel &radModel,
                    vector<double> &A, vector<double> &B, vector<double> &C)
    {
        const int &nZones = radModel.params.nZones;
        const int &nQuad = radModel.params.nQuad;
        const int halfQuad = int(0.5 * nQuad);
        const vector<double> &tau = radModel.tau;
        const vector<double> &quadMu = radModel.quadMu;
        const vector<double> &quadW = radModel.quadW;
        double quadSum;

        // Yes they are nZones even though I only use nZones - 1 or - 2 #YOLO
        vector<vector<double>> Dtau(nZones, vector<double>(halfQuad));
        vector<vector<double>> expDtau(nZones, vector<double>(halfQuad));

        vector<vector<double>> e0(nZones, vector<double>(halfQuad));
        vector<vector<double>> e1(nZones, vector<double>(halfQuad));
        vector<vector<double>> e2(nZones, vector<double>(halfQuad));

        vector<vector<double>> alphaP(nZones, vector<double>(halfQuad));
        vector<vector<double>> betaP(nZones, vector<double>(halfQuad));
        vector<vector<double>> gammaP(nZones, vector<double>(halfQuad));

        vector<vector<double>> alphaM(nZones, vector<double>(halfQuad));
        vector<vector<double>> betaM(nZones, vector<double>(halfQuad));
        vector<vector<double>> gammaM(nZones, vector<double>(halfQuad));

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
        for (int i=1; i < nZones; i++)
        {
            for (int j=0; j < halfQuad; j++)
            {
                alphaM[i][j] = e0[i][j] +
                               (e2[i][j] - (Dtau[i][j] + 2 * Dtau[i - 1][j]) * e1[i][j]) /
                               (Dtau[i - 1][j] * (Dtau[i][j] + Dtau[i - 1][j]));
                betaM[i][j] = ((Dtau[i][j] + Dtau[i - 1][j]) * e1[i][j] - e2[i][j]) /
                              ((Dtau[i][j] * Dtau[i - 1][j]));
                gammaM[i][j] = (e2[i][j] - Dtau[i - 1][j] * e1[i][j]) /
                               (Dtau[i][j] * (Dtau[i][j] + Dtau[i - 1][j]));
            }
        }

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

        // Integrate mu first (Quadrature)

        // C integral
        for (int i=1; i < nZones - 1; i++)
        {
            quadSum = zero;
            for (int j=0; j < halfQuad; j++)
            {
                quadSum += quadW[j] * (gammaM[i - 1][j] +
                                       gammaP[i - 1][j] + expDtau[i - 1][j] *
                                                          (alphaP[i + 1][j] * expDtau[i][j] + betaP[i][j]));
            }

            C[i] = 0.5 * quadSum;
        }

        // B integral
        for (int i=1; i < nZones - 1; i++)
        {
            quadSum = zero;
            for (int j=0; j < halfQuad; j++)
            {
                quadSum += quadW[j] * (betaM[i][j] + gammaM[i - 1][j] * expDtau[i - 1][j] +
                                       betaP[i][j] + alphaP[i + 1][j] * expDtau[i][j]);
            }

            B[i] = 0.5 * quadSum;
        }

        // A integral
        for (int i=1; i < nZones - 1; i++)
        {
            quadSum = zero;
            for (int j=0; j < halfQuad; j++)
            {
                quadSum += quadW[j] * (alphaM[i + 1][j] + expDtau[i][j] *
                                                          (gammaM[i - 1][j] * expDtau[i - 1][j] + betaM[i][j]) +
                                       alphaP[i + 1][j]);
            }

            A[i] = 0.5 * quadSum;
        }

        // 4 unknown matrix elements at i = 1 and i = N
        quadSum = zero;
        for (int j=0; j < halfQuad; j++)
        {
            quadSum += quadW[j] * (betaP[nZones - 1][j] * expDtau[nZones - 2][j] + gammaP[nZones - 2][j]);
        }
        C[nZones - 1] = quadSum;

        // quadSum = zero;
        // for (int j=0; j < halfQuad; j++)
        // {
        //     quadSum += quadW[j] * (e1[nZones - 1][j] / Dtau[nZones - 2][j] +
        //                            (Dtau[nZones - 1][j] - 1 + expDtau[nZones - 1][j]) / Dtau[nZones - 1][j]);
        // }
        B[nZones - 1] = 1.0;

        quadSum = zero;
        for (int j=0; j < halfQuad; j++)
        {
            // Remember i^-(0) = 0 so this is just i^+(0) in integral
            quadSum += quadW[j] * (e1[1][j] / Dtau[0][j]);
        }
        B[0] = quadSum;

        quadSum = zero;
        for (int j=0; j < halfQuad; j++)
        {
            // Remember i^-(0) = 0 so i^-_{1 + 1} = alphaM_{1 + 1}
            quadSum += quadW[j] * (alphaM[1][j]);
        }
        A[0] = quadSum;

        // Finally integrate across wavelengths
        for (int i=1; i < nZones - 1; i++)
        {

        }
    }
}