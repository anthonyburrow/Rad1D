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
        const double &nZones = nuModel.params.nZones;
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
        const double &nZones = nuModel.params.nZones;
        const double &eps = nuModel.params.eps;
        const vector<vector<double>> &lambda = nuModel.lambda;
        const vector<double> &lambdaA = nuModel.lambdaA;
        const vector<double> &lambdaB = nuModel.lambdaB;
        const vector<double> &lambdaC = nuModel.lambdaC;

        vector<double> x(nZones, 0.0);
        vector<double> y(nZones, 0.0);
        vector<double> JFormal(nZones, 0.0);

        for (int i = 0; i < nZones; i++)
        {
            for (int j = 0; j < nZones; j++)
            {
                JFormal[i] += lambda[i][j] * nuModel.S[j];
            }
            x[i] = JFormal[i] - nuModel.J[i];
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
}