#include <vector>
#include <cmath>
#include <iostream>

#include "lambdaIterate.hpp"
#include "NuModel.hpp"
#include "constants.hpp"
#include "blackbody.hpp"
#include "util.hpp"

using namespace std;

namespace myLib
{
    void calcJ(NuModel &nuModel)
    {
        const double &nZones = nuModel.params.nZones;
        const vector<vector<double>> &lambda = nuModel.lambda;
        double jTerm;

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
            nuModel.S[i] = eps * nuModel.B[i] + (1.0 - eps) * nuModel.J[i];
        }
    }


    //SMS moronic attempts at coding
    void ALIcalcJ(NuModel &nuModel)
    {
        const double &nZones = nuModel.params.nZones;
        const double &eps = nuModel.params.eps;
        const vector<vector<double>> &lambda = nuModel.lambda;
        double jTerm;
        
        std::vector<double> Jold = nuModel.J;
        std::vector<double> x(nZones,0.0);
        std::vector<double> y(nZones,0.0);
        std::vector<double> formal_J(nZones,0.0);

        for (int i=0; i < nZones; i++)
        {
            jTerm = 0.0;
            for (int j=0; j < nZones; j++)
            {
                jTerm += lambda[i][j] * nuModel.S[j];
            }
            formal_J[i] = jTerm;
            x[i] = formal_J[i] - Jold[i];
        }
        
        TridiagonalSoln(y, eps, lambda, x);
        
        for (int i=0; i < nZones; i++)
        {
            nuModel.J[i] = y[i] + Jold[i];
        }        
    }

    //SMS moronic attempts at coding
    void ALIcalcS(NuModel &nuModel)
    {
        const int &nZones = nuModel.params.nZones;
        const double &eps = nuModel.params.eps;

        const vector<vector<double>> &lambda = nuModel.lambda;
        std::vector<double> Sold = nuModel.S;

        double sTerm;

        for (int i=0; i < nZones; i++)
        {
            sTerm = 0.0;
            for (int j=0; j < nZones; j++)
            {
                sTerm += (1.0 - eps) * lambda[i][j] * Sold[j];
            }

            nuModel.S[i] = eps * nuModel.B[i] +  sTerm;
        }

    }

    void lambdaIteration(NuModel &nuModel)
    {
        calcJ(nuModel);
        calcS(nuModel);
    }

    //SMS moronic attempts at coding
    void ALI(NuModel &nuModel)
    {
        ALIcalcJ(nuModel);
        // ALIcalcS(nuModel);
        calcS(nuModel);
    }
}