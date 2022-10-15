#include <iostream>

#include "NuModel.hpp"
#include "constants.hpp"
#include "process.hpp"
#include "blackbody.hpp"

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
        // I(params.nZones, vector<double>(params.nQuad, zero)),
        J(params.nZones, zero)
    {
        setBoundary();
    }

    void NuModel::setBoundary()
    {
        const int &nZones = params.nZones;
        // const int &nQuad = params.nQuad;
        // const int halfQuad = int(0.5 * nQuad);
        double bb;

        // Initialize { I+(tauMax) = B(T) ; I-(0) = 0 }
        // for (int i = 0; i < halfQuad; i++)
        // {
        //     I[0][halfQuad + i] = zero;
        //     I[nZones - 1][i] = planck(lam, radModel.T[nZones - 1]);
        // }

        // Initialize S(tau) = B(T(tau))
        for (int i = 0; i < nZones; i++)
        {
            bb = planck(lam, radModel.T[i]);
            B[i] = bb;
            S[i] = bb;
        }
    }

    const double NuModel::getFlux()
    {
        return calcFlux(*this);
    }
}
