#include <vector>
#include <cmath>
#include <iostream>

#include "lambdaIterate.hpp"
#include "NuModel.hpp"
#include "constants.hpp"
#include "blackbody.hpp"
#include "util.hpp"

using namespace std;

namespace myLib
{
    void calcJ(NuModel &nuModel)
    {
        const int &nZones = nuModel.params.nZones;
        const vector<vector<double>> &lambda = nuModel.lambda;
        double jTerm;

        for (int i = 0; i < nZones; i++)
        {
            jTerm = 0.0;
            for (int j = 0; j < nZones; j++)
            {
                jTerm += lambda[i][j] * nuModel.S[j];
            }

            nuModel.J[i] = jTerm;
        }
    }

    void ALIcalcJ(NuModel &nuModel)
    {
        const int &nZones = nuModel.params.nZones;
        const double &eps = nuModel.params.eps;
        const vector<vector<double>> &lambda = nuModel.lambda;
        const vector<double> &lambdaA = nuModel.lambdaA;
        const vector<double> &lambdaB = nuModel.lambdaB;
        const vector<double> &lambdaC = nuModel.lambdaC;

        vector<double> x(nZones, 0.0);
        vector<double> y(nZones, 0.0);
        double JFormal;

        for (int i = 0; i < nZones; i++)
        {
            JFormal = 0.0;
            for (int j = 0; j < nZones; j++)
            {
                JFormal += lambda[i][j] * nuModel.S[j];
            }
            x[i] = JFormal - nuModel.J[i];
        }

        TridiagonalSoln(y, lambdaA, lambdaB, lambdaC, x);

        for (int i = 0; i < nZones; i++)
        {
            nuModel.J[i] += y[i];
        }        
    }

    void calcS(NuModel &nuModel)
    {
        const int &nZones = nuModel.params.nZones;
        const double &eps = nuModel.params.eps;

        for (int i = 0; i < nZones; i++)
        {
            nuModel.S[i] = eps * nuModel.B[i] + (1.0 - eps) * nuModel.J[i];
        }
    }

    void lambdaIteration(NuModel &nuModel)
    {
        calcJ(nuModel);
        calcS(nuModel);
    }

    void ALI(NuModel &nuModel)
    {
        ALIcalcJ(nuModel);
        calcS(nuModel);
    }

    void NgIteration(NuModel &nuModel, const vector<double> &S3,
                     const vector<double> &S2, const vector<double> &S1,
                     const vector<double> &S0)
    {
        const int &nZones = nuModel.params.nZones;
        const double &eps = nuModel.params.eps;
        vector<double> &Snew = nuModel.S;
        vector<double> &Jnew = nuModel.J;

        double Q1, Q2, Q3;

        double A1 = 0.0, A2 = 0.0;
        double &B1 = A2;
        double B2 = 0.0;
        double C1 = 0.0, C2 = 0.0;

        double a, b;

        // Get Snew using S3 - S0 here
        for (int i = 0; i < nZones; i++)
        {
            Q1 = S0[i] - 2.0 * S1[i] + S2[i];
            Q2 = S0[i] - S1[i] - S2[i] + S3[i];
            Q3 = S0[i] - S1[i];

            A1 += Q1 * Q1;
            A2 += Q2 * Q1;
            B2 += Q2 * Q2;
            C1 += Q1 * Q3;
            C2 += Q2 * Q3;
        }

        a = (C1 * B2 - C2 * B1) / (A1 * B2 - A2 * B1);
        b = (C2 * A1 - C1 * A2) / (A1 * B2 - A2 * B1);

        for (int i = 0; i < nZones; i++)
        {
            Snew[i] = (1.0 - a - b) * S0[i] + a * S1[i] + b * S2[i];
            Jnew[i] = (Snew[i] - eps * nuModel.B[i]) / (1.0 - eps);
        }
    }
}