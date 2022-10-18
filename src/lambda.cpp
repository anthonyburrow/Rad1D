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

        // Yes they are nZones even though I sometimes use nZones - 1 or - 2 #YOLO
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

    void calcLambda(const RadModel &radModel, vector<vector<double>> &lambda)
    {
        const int &nZones = radModel.params.nZones;
        const int &nQuad = radModel.params.nQuad;
        const int halfQuad = int(0.5 * nQuad);
        const vector<double> &tau = radModel.tau;
        const vector<double> &quadMu = radModel.quadMu;
        const vector<double> &quadW = radModel.quadW;
        double quadSum;
        double Ip;
        double Im;
        double ImPrev;
        double IpPrev;

        // Yes they are nZones even though I sometimes use nZones - 1 or - 2 #YOLO
        vector<vector<double>> Dtau(nZones, vector<double>(halfQuad, zero));
        vector<vector<double>> expDtau(nZones, vector<double>(halfQuad, zero));

        vector<vector<double>> e0(nZones, vector<double>(halfQuad, zero));
        vector<vector<double>> e1(nZones, vector<double>(halfQuad, zero));
        vector<vector<double>> e2(nZones, vector<double>(halfQuad, zero));

        vector<vector<double>> alphaP(nZones, vector<double>(halfQuad, zero));
        vector<vector<double>> betaP(nZones, vector<double>(halfQuad, zero));
        vector<vector<double>> gammaP(nZones, vector<double>(halfQuad, zero));

        vector<vector<double>> alphaM(nZones, vector<double>(halfQuad, zero));
        vector<vector<double>> betaM(nZones, vector<double>(halfQuad, zero));
        vector<vector<double>> gammaM(nZones, vector<double>(halfQuad, zero));

        // Calc delta tau's
        for (int i=1; i < nZones; i++)
        {
            for (int j=0; j < halfQuad; j++)
            {
                Dtau[i - 1][j] = (tau[i] - tau[i - 1]) / quadMu[j];
                expDtau[i - 1][j] = exp(-Dtau[i - 1][j]);
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

        for (int j=0; j < halfQuad; j++)
        {
            // Lol i have no idea for these...
            // I just used the linear interpolation values for these
            alphaM[nZones - 1][j] = e0[nZones - 1][j] - e1[nZones - 1][j] / Dtau[nZones - 2][j];   
            alphaP[nZones - 1][j] = zero;

            betaP[0][j] = e1[1][j] / Dtau[0][j];
            betaM[nZones - 1][j] = e1[nZones - 1][j] / Dtau[nZones - 2][j];

            gammaM[0][j] = zero;
            gammaP[0][j] = e0[1][j] - e1[1][j] / Dtau[0][j];
        }

        // Integrate mu first (Quadrature) to fill Lambda

        // Lambda columns 2 to (N - 1)
        for (int i=1; i < nZones - 1; i++)
        {
            // C (i - 1, i)
            quadSum = zero;
            for (int j=0; j < halfQuad; j++)
            {
                Im = gammaM[i - 1][j];
                Ip = (alphaP[i + 1][j] * expDtau[i][j] + betaP[i][j]) * expDtau[i - 1][j] +
                     gammaP[i - 1][j];
                quadSum += quadW[j] * (Im + Ip);
            }
            lambda[i - 1][i] = 0.5 * quadSum;

            // B (i, i)
            quadSum = zero;
            for (int j=0; j < halfQuad; j++)
            {
                Im = gammaM[i - 1][j] * expDtau[i - 1][j] + betaM[i][j];
                Ip = alphaP[i + 1][j] * expDtau[i][j] + betaP[i][j];
                quadSum += quadW[j] * (Im + Ip);
            }
            lambda[i][i] = 0.5 * quadSum;

            // A (i + 1, i)
            quadSum = zero;
            for (int j=0; j < halfQuad; j++)
            {
                Im = (gammaM[i - 1][j] * expDtau[i - 1][j] + betaM[i][j]) * expDtau[i][j] + alphaM[i + 1][j];
                Ip = alphaP[i + 1][j];
                quadSum += quadW[j] * (Im + Ip);
            }
            lambda[i + 1][i] = 0.5 * quadSum;

            // I^- from row (i + 2) to N; I^+ = 0 in this section
            for (int j=0; j < halfQuad; j++)
            {
                // This is the A I^- calculation
                ImPrev = (gammaM[i - 1][j] * expDtau[i - 1][j] + betaM[i][j]) * expDtau[i][j] + alphaM[i + 1][j];
                for (int k=i + 2; k < nZones; k++)
                {
                    Im = ImPrev * expDtau[k - 1][j];
                    lambda[k][i] += 0.5 * quadW[j] * Im;
                    ImPrev = Im;
                }
            }

            // I^+ from row (i - 2) to 1; I^- = 0 in this section
            for (int j=0; j < halfQuad; j++)
            {
                // This is the C I^+ calculation
                IpPrev = (alphaP[i + 1][j] * expDtau[i][j] + betaP[i][j]) * expDtau[i - 1][j] +
                         gammaP[i - 1][j];
                for (int k=i - 1; k > -1; k--)
                {
                    Ip = IpPrev * expDtau[k][j];
                    lambda[k][i] += 0.5 * quadW[j] * Ip;
                    IpPrev = Ip;
                }
            }
        }

        // Column 1
        quadSum = zero;
        for (int j=0; j < halfQuad; j++)
        {
            // i^-(0) = 0 so this is just i^+(0) in integral
            quadSum += quadW[j] * betaP[0][j];
        }
        lambda[0][0] = 0.5 * quadSum;

        quadSum = zero;
        for (int j=0; j < halfQuad; j++)
        {
            // i^-(0) = 0 so i^-_{1 + 1} = alphaM_{1 + 1}
            quadSum += quadW[j] * alphaM[1][j];
        }
        lambda[1][0] = 0.5 * quadSum;

        // I^- from row 3 to N; I^+ = 0 in this section
        for (int j=0; j < halfQuad; j++)
        {
            // This is the lambda[1][0] I^- calculation above
            ImPrev = alphaM[1][j];
            for (int k=2; k < nZones; k++)
            {
                Im = ImPrev * expDtau[k - 1][j];
                lambda[k][0] += 0.5 * quadW[j] * Im;
                ImPrev = Im;
            }
        }

        // Column N
        lambda[nZones - 1][nZones - 1] = 1.0;   // Not sure about this one

        quadSum = zero;
        for (int j=0; j < halfQuad; j++)
        {
            Ip = 1.0 * expDtau[nZones - 2][j] + gammaP[nZones - 2][j];
            quadSum += quadW[j] * Ip;
        }
        lambda[nZones - 2][nZones - 1] = 0.5 * quadSum;

        // I^+ from row (N - 2) to 1; I^- = 0 in this section
        for (int j=0; j < halfQuad; j++)
        {
            // This is the C I^+ calculation
            IpPrev = 1.0 * expDtau[nZones - 2][j] + gammaP[nZones - 2][j];
            for (int k=nZones - 3; k > -1; k--)
            {
                Ip = IpPrev * expDtau[k][j];
                lambda[k][nZones - 1] += 0.5 * quadW[j] * Ip;
                IpPrev = Ip;
            }
        }
    }
}