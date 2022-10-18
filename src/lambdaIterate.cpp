#include <vector>
#include <cmath>
#include <iostream>

#include "lambdaIterate.hpp"
#include "NuModel.hpp"
#include "constants.hpp"
#include "blackbody.hpp"

using namespace std;

namespace myLib
{
    void calcJ(NuModel &nuModel)
    {
        const double &nZones = nuModel.params.nZones;
        // const vector<double> &lambdaA = nuModel.radModel.lambdaA;
        // const vector<double> &lambdaB = nuModel.radModel.lambdaB;
        // const vector<double> &lambdaC = nuModel.radModel.lambdaC;
        const vector<vector<double>> &lambda = nuModel.radModel.lambda;

        nuModel.J[nZones - 1] = nuModel.S[nZones - 1];

        for (int i=0; i < nZones - 1; i++)
        {
            for (int j=0; j < nZones; j++)
            {
                nuModel.J[i] += lambda[i][j] * nuModel.S[j];
            }
        }
    }

    void calcS(NuModel &nuModel)
    {
        const int &nZones = nuModel.params.nZones;
        const double &eps = nuModel.params.eps;

        for (int i=0; i < nZones; i++)
        {
            nuModel.S[i] = eps * nuModel.B[i] + (1.0 - eps) * nuModel.J[i];
        }
    }

    void lambdaIteration(NuModel &nuModel)
    {
        calcJ(nuModel);
        calcS(nuModel);
    }
}