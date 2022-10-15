#include <vector>
#include <cmath>
#include <iostream>

#include "lambdaIterate.hpp"
#include "NuModel.hpp"
#include "process.hpp"
#include "constants.hpp"
#include "blackbody.hpp"

using namespace std;

namespace myLib
{
    // void calcJ(NuModel &nuModel)
    // {
    //     const int &nZones = nuModel.params.nZones;
    //     const int &nQuad = nuModel.params.nQuad;
    //     double J;

    //     for (int i=0; i < nZones; i++)
    //     {
    //         J = 0.0;
    //         for (int j=0; j < nQuad; j++)
    //         {
    //             // J_i = I_ij W_j
    //             J += nuModel.I[i][j] * nuModel.radModel.quadW[j];
    //         }

    //         nuModel.J[i] = 0.5 * J;
    //     }
    // }

    void calcS(NuModel &nuModel)
    {
        const int &nZones = nuModel.params.nZones;
        const double &eps = nuModel.params.eps;

        for (int i=0; i < nZones; i++)
        {
            nuModel.S[i] = eps * planck(nuModel.lam, nuModel.radModel.T[i]) +
                           (1.0 - eps) * nuModel.J[i];
        }
    }

    void lambdaIteration(NuModel &nuModel)
    {
        // calcJ(nuModel);
        // calcS(nuModel);
    }
}