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
        const vector<double> &lambdaA = nuModel.radModel.lambdaA;
        const vector<double> &lambdaB = nuModel.radModel.lambdaB;
        const vector<double> &lambdaC = nuModel.radModel.lambdaC;

        nuModel.J[0] = lambdaB[0] * nuModel.S[0] + lambdaC[1] * nuModel.S[1];
        // nuModel.J[nZones - 1] = lambdaA[nZones - 2] * nuModel.S[nZones - 2] +
        //                         lambdaB[nZones - 1] * nuModel.S[nZones - 1];
        nuModel.J[nZones - 1] = nuModel.S[nZones - 1];

        for (int i=1; i < nZones - 1; i++)
        {
            nuModel.J[i] = lambdaA[i - 1] * nuModel.S[i - 1] +
                           lambdaB[i] * nuModel.S[i] +
                           lambdaC[i + 1] * nuModel.S[i + 1];
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