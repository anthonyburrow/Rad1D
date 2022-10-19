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
        double quadSum;
        double Ip;
        double Im;
        double deltaIp;
        double deltaIm;
        double ImPrev;
        double IpPrev;

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

                // Can't be < 0 because it is an intensity value
                gammaM[i][j] = max(gammaM[i][j], zero);
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

                // Can't be < 0 because it is an intensity value
                alphaP[i][j] = max(alphaP[i][j], zero);
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

            // I^+(tau_max) normalized
            betaP[nZones - 1][j] = 1.0;

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
                deltaIm = max(gammaM[i - 1][j], zero);
                deltaIp = max(gammaP[i - 1][j], zero);

                Im = deltaIm;
                Ip = (alphaP[i + 1][j] * expDtau[i][j] + betaP[i][j]) * expDtau[i - 1][j] +
                     deltaIp;
                quadSum += quadW[j] * (Im + Ip);
            }
            lambda[i - 1][i] = 0.5 * quadSum;

            // B (i, i)
            quadSum = zero;
            for (int j=0; j < halfQuad; j++)
            {
                deltaIm = max(betaM[i][j], zero);
                deltaIp = max(betaP[i][j], zero);

                Im = gammaM[i - 1][j] * expDtau[i - 1][j] + deltaIm;
                Ip = alphaP[i + 1][j] * expDtau[i][j] + deltaIp;
                quadSum += quadW[j] * (Im + Ip);
            }
            lambda[i][i] = 0.5 * quadSum;

            // A (i + 1, i)
            quadSum = zero;
            for (int j=0; j < halfQuad; j++)
            {
                deltaIm = max(alphaM[i + 1][j], zero);
                deltaIp = max(alphaP[i + 1][j], zero);

                Im = (gammaM[i - 1][j] * expDtau[i - 1][j] + betaM[i][j]) * expDtau[i][j] + deltaIm;
                Ip = deltaIp;
                quadSum += quadW[j] * (Im + Ip);
            }
            lambda[i + 1][i] = 0.5 * quadSum;

            // I^- from row (i + 2) to N; I^+ = 0 in this section
            for (int j=0; j < halfQuad; j++)
            {
                // This is the A I^- calculation
                deltaIm = max(alphaM[i + 1][j], zero);
                ImPrev = (gammaM[i - 1][j] * expDtau[i - 1][j] + betaM[i][j]) * expDtau[i][j] +
                         deltaIm;
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
                deltaIp = max(gammaP[i - 1][j], zero);
                IpPrev = (alphaP[i + 1][j] * expDtau[i][j] + betaP[i][j]) * expDtau[i - 1][j] +
                         deltaIp;
                for (int k=i - 2; k > -1; k--)
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
            deltaIp = max(betaP[0][j], zero);
            quadSum += quadW[j] * deltaIp;
        }
        lambda[0][0] = 0.5 * quadSum;

        quadSum = zero;
        for (int j=0; j < halfQuad; j++)
        {
            // i^-(0) = 0 so i^-_{1 + 1} = alphaM_{1 + 1}
            deltaIm = max(alphaM[1][j], zero);
            quadSum += quadW[j] * deltaIm;
        }
        lambda[1][0] = 0.5 * quadSum;

        // I^- from row 3 to N; I^+ = 0 in this section
        for (int j=0; j < halfQuad; j++)
        {
            // This is the lambda[1][0] I^- calculation above
            deltaIm = max(alphaM[1][j], zero);
            ImPrev = deltaIm;
            for (int k=2; k < nZones; k++)
            {
                Im = ImPrev * expDtau[k - 1][j];
                lambda[k][0] += 0.5 * quadW[j] * Im;
                ImPrev = Im;
            }
        }

        // Column N
        quadSum = zero;
        for (int j=0; j < halfQuad; j++)
        {
            deltaIm = max(betaM[nZones - 1][j], zero);
            Im = deltaIm;
            Ip = betaP[nZones - 1][j];
            quadSum += quadW[j] * (Im + Ip);
        }
        lambda[nZones - 1][nZones - 1] = 0.5 * quadSum;

        quadSum = zero;
        for (int j=0; j < halfQuad; j++)
        {
            // Here assume I^+(tau_max) = 1 (normalized)?
            deltaIp = max(gammaP[nZones - 2][j], zero);
            Ip = betaP[nZones - 1][j] * expDtau[nZones - 2][j] + deltaIp;
            quadSum += quadW[j] * Ip;
        }
        lambda[nZones - 2][nZones - 1] = 0.5 * quadSum;

        // I^+ from row (N - 2) to 1; I^- = 0 in this section
        for (int j=0; j < halfQuad; j++)
        {
            // This is the C I^+ calculation
            deltaIp = max(gammaP[nZones - 2][j], zero);
            IpPrev = betaP[nZones - 1][j] * expDtau[nZones - 2][j] + deltaIp;
            for (int k=nZones - 3; k > -1; k--)
            {
                Ip = IpPrev * expDtau[k][j];
                lambda[k][nZones - 1] += 0.5 * quadW[j] * Ip;
                IpPrev = Ip;
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