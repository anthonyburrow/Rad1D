#include <iostream>

#include "NuModel.hpp"
#include "constants.hpp"
#include "blackbody.hpp"
#include "lambdaIterate.hpp"

using namespace std;

namespace myLib
{
    NuModel::NuModel(const double &lam, const RadModel &radModel)
      : lam(lam),
        radModel(radModel),
        params(radModel.params),
        // Allocate vectors
        B(params.nZones, zero),
        S(params.nZones, zero),
        J(params.nZones, zero)
    {
        setBoundary();
    }

    void NuModel::setBoundary()
    {
        const int &nZones = params.nZones;
        double bb;

        // Initialize S(tau) = B(T(tau))
        for (int i = 0; i < nZones; i++)
        {
            bb = planck(lam, radModel.T[i]);
            B[i] = bb;
            S[i] = bb;
        }
    }

    void NuModel::iterate()
    {
        lambdaIteration(*this);
    }

    const double NuModel::calcFlux()
    {
        const int &maxIter = params.maxIter;
        // const int &nQuad = nuModel.params.nQuad;
        double F = 0.0;

        // Converge S & J
        for (int i=0; i < maxIter; i++)
        {
            iterate();

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
