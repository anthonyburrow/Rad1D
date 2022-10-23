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
        const vector<vector<double>> &lambda = nuModel.radModel.lambda;
        double jTerm;

        // nuModel.J[nZones - 1] = nuModel.S[nZones - 1];

        // for (int i=0; i < nZones - 1; i++)
        for (int i=0; i < nZones; i++)
        {
            jTerm = 0.0;
            for (int j=0; j < nZones; j++)
            {
                jTerm += lambda[i][j] * nuModel.S[j];
            }

            nuModel.J[i] = jTerm;
        }
    }

    void calcS(NuModel &nuModel)
    {
        const int &nZones = nuModel.params.nZones;
        const double &eps = nuModel.params.eps;

        for (int i=0; i < nZones; i++)
        {
            nuModel.S[i] = eps + (1.0 - eps) * nuModel.J[i];
        }
    }

    void lambdaIteration(NuModel &nuModel)
    {
        calcJ(nuModel);
        calcS(nuModel);
    }
}