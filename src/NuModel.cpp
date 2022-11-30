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
        // Physical quantites
        tau(params.nZones, zero),
        B(params.nZones, zero),
        S(params.nZones, zero),
        J(params.nZones, zero),
        // Lambda iteration vectors
        lambda(params.nZones, vector<double>(params.nZones, 0.0)),
        lambdaA(params.nZones - 1, zero),
        lambdaB(params.nZones, zero),
        lambdaC(params.nZones - 1, zero)
    {
        setInitialCond();
        calcTau();
        initLambda();
    }

    void NuModel::setInitialCond()
    {
        const int &nZones = params.nZones;

        // Initialize S(tau) = B(T(tau))
        for (int i = 0; i < nZones; i++)
        {
            B[i] = 1.0;
            // B[i] = planck(lam, radModel.T[i]);
            S[i] = B[i];
        }
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

    void NuModel::initLambda()
    {
        const int &nZones = params.nZones;
        const double &eps = params.eps;

        calcLambda(*this);

        if (!params.accelerated) { return; }

        for (int i = 0; i < nZones; i++)
        {
            lambdaB[i] = 1.0 - (1.0 - eps) * lambda[i][i];
        }
        for(int i = 0; i < nZones - 1; i++)
        {
            lambdaA[i] = -(1.0 - eps) * lambda[i + 1][i];
            lambdaC[i] = -(1.0 - eps) * lambda[i][i + 1];
        }
    }

    void NuModel::iterate(const bool accelerate)
    {
        if (accelerate && params.accelerated) { ALI(*this); }
        else { lambdaIteration(*this); }
    }

    double NuModel::calcFlux()
    {
        const int &maxIter = params.maxIter;
        const double &epsConverge = params.epsConverge;
        double flux;
        double prevJ = 1;

        // Do Ng Acceleration here

        // Converge S & J
        iterate(false);

        for (int i = 2; i < maxIter; i++)
        {
            iterate();

            // Check/break for convergence at the surface
            if (2.0 * abs(J[0] - prevJ) / (J[0] + prevJ) < epsConverge) { break; }
            prevJ = J[0];
        }

        // Calc F based on converged S & J
        flux = calcF0();

        return flux;
    }

    double NuModel::calcF0()
    {
        const int &nZones = params.nZones;
        const int &nQuad = params.nQuad;
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

        for (int j = 0; j < nQuad; j++)
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

        return F0;
    }
}
