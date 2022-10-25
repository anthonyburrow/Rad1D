#include <iostream>
#include <cmath>

#include "NuModel.hpp"
#include "constants.hpp"
#include "blackbody.hpp"
#include "lambdaIterate.hpp"
#include "lambda.hpp"

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
        const double bbScale = 2.0 * hc * c * 1e8;

        // Initialize S(tau) = B(T(tau)) = 1
        for (int i = 0; i < nZones; i++)
        {
            B[i] = bbScale * planck(lam, radModel.T[i]);
            S[i] = 1.0;
        }
    }

    void NuModel::iterate()
    {
        lambdaIteration(*this);
    }

    double NuModel::calcFlux()
    {
        const int &maxIter = params.maxIter;
        double flux;

        // Converge S & J
        for (int i=0; i < maxIter; i++)
        {
            iterate();

            // Check/break for convergence?
        }

        // Calc F based on converged S & J
        flux = calcF0();
        flux *= B[0];

        return flux;
    }

    double NuModel::calcF0()
    {
        const int &nZones = params.maxIter;
        const int halfQuad = int(0.5 * params.nQuad);
        const vector<double> &tau = radModel.tau;
        const vector<double> &quadMu = radModel.quadMu;
        const vector<double> &quadW = radModel.quadW;

        double F0 = zero;
        double Ip0;

        vector<double> Dtau(nZones, zero);
        vector<double> expDtau(nZones, zero);

        vector<double> e0(nZones, zero);
        vector<double> e1(nZones, zero);
        vector<double> e2(nZones, zero);

        vector<vector<double>> Ip(nZones, vector<double>(nZones, zero));

        for (int j = 0; j < halfQuad; j++)
        {
            // Setup values needed for both I+ and I-
            setupFormalSoln(Dtau, expDtau, e0, e1, e2,
                            quadMu[j], tau, nZones);

            // Solve for I+ matrix (I^-(0) = 0)
            calcIpMatrix(Ip, Dtau, expDtau, e0, e1, e2, nZones);

            Ip0 = zero;
            for (int i = 0; i < nZones; i++)
            {
                Ip0 += Ip[0][i] * S[i];
            }
            F0 += quadW[j] * quadMu[j] * Ip0;
        }

        // F(0) = 4 pi H(0) = 2 pi mu_j I^+_j(0) Wj
        F0 *= 2.0 * pi;

        // cout << F0 << endl;

        return F0;
    }
}
