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

    void NgIteration(NuModel &nuModel, vector<double> &S3, vector<double> &S2,
                     vector<double> &S1, vector<double> &S0)
    {
        const int &nZones = nuModel.params.nZones;
        vector<double> &Snew = nuModel.S;

        // Get Snew using S3 - S0 here
        for( int i = 0; i < nZones; i++)
        {

            double Q1 = S0[i] - 2.*S1[i] + S2[i];
            double Q2 = S0[i] - S1[i] - S2[i] + S3[i];
            double Q3 = S0[i] - S1[i];

            double A1 = Q1 * Q1;
            double A2 = Q2 * Q1;

            double B1 = Q1 * Q2;
            double B2 = Q2 * Q2;

            double C1 = Q1 * Q3;
            double C2 = Q2 * Q3;

            double a = (C1*B2 - C2*B2)/(A1*B2 - A2*B1);
            double b = (C2*A1 - C1*A2)/(A1*B2 - A2*B1);

            Snew[i] = (1. - a - b)*S0[i]+a*S1[i]+b*S2[i];

        }

        nuModel.S = Snew;
    }
}