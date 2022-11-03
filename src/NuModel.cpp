#include <iostream>
#include <cmath>

#include "NuModel.hpp"
#include "constants.hpp"
#include "blackbody.hpp"
#include "lambdaIterate.hpp"
#include "lambda.hpp"
#include "lineProfiles.hpp"
#include "io.hpp"

using namespace std;

namespace myLib
{
    NuModel::NuModel(const double &lam, const RadModel &radModel)
      : lam(lam),
        radModel(radModel),
        params(radModel.params),
        // Allocate vectors
        tau(params.nZones, zero),
        B(params.nZones, zero),
        S(params.nZones, zero),
        J(params.nZones, zero),
        lambda(params.nZones, vector<double>(params.nZones, 0.0))
    {
        setInitialCond();
        calcTau();
        calcLambda(*this, lambda);
    }

    void NuModel::calcTau()
    {
        const int &nZones = params.nZones;
        double frac;

        for (int i = 0; i < nZones; i++)
        {
            frac = zero;
            for (feature line : radModel.lineList)
            {
                frac += line.tauRatio * gaussianProfile(lam, radModel.T[i], line);
            }
            tau[i] = radModel.tauCont[i] * (1.0 + frac);
        }
    }

    void NuModel::setInitialCond()
    {
        const int &nZones = params.nZones;
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
        const double &epsConverge = params.epsConverge;
        double flux;

        double prevJ = 1;

        // Do Ng Acceleration here

        // Converge S & J
        for (int i=0; i < maxIter; i++)
        {
            iterate();

            // Check/break for convergence at the surface
            if (2.0 * abs(J[0] - prevJ) / (J[0] + prevJ) < epsConverge) { break; }
            prevJ = J[0];
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
                // Check for NaN (not sure why there are NaNs in this matrix...)
                if (Ip[0][i] != Ip[0][i]) { continue; }

                Ip0 += Ip[0][i] * S[i];
            }
            F0 += quadW[j] * quadMu[j] * Ip0;
        }

        // F(0) = 4 pi H(0) = 2 pi mu_j I^+_j(0) Wj
        F0 *= 2.0 * pi;

        return F0;
    }
}
