#include <vector>
#include <cmath>
#include <iostream>

#include "NuModel.hpp"
#include "process.hpp"
#include "constants.hpp"
#include "lambdaIterate.hpp"

using namespace std;

namespace myLib
{
    double calcFlux(NuModel &nuModel)
    {
        const int &maxIter = nuModel.params.maxIter;
        // const int &nQuad = nuModel.params.nQuad;
        double F = 0.0;

        // Converge S & J
        for (int i=0; i < maxIter; i++)
        {
            lambdaIteration(nuModel);

            // Check/break for convergence?
        }

        // Calc F based on converged S/J/I
        // F(0) = 4 pi H(0) = 2 pi mu_j I_0j Wj
        // for (int i=0; i < nQuad; i++)
        // {
        //     F += nuModel.radModel.quadMu[i] * nuModel.I[0][i] * nuModel.radModel.quadW[i];
        // }

        // F *= 2.0 * pi;

        return F;
    }
}