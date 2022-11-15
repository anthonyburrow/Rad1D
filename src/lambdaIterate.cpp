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
        std::vector<double> a(nZones-1,0.0);
        std::vector<double> b(nZones,0.0);
        std::vector<double> c(nZones-1,0.0);
        std::vector<double> x(nZones,0.0);
        std::vector<double> y(nZones,0.0);
        std::vector<double> formal_J(nZones,0.0);

        for (int i=0; i < nZones; i++)
        {
            jTerm = 0.0;
            for (int j=0; j < nZones; j++)
            {
                jTerm += lambda[i][j] * nuModel.S[j];
                if(i==j)
                {
                    b[i]=1.0-(1.0-eps)*lambda[i][j];
                }
                if(i==j+1)
                {
                    a[i]=1.0-(1.0-eps)*lambda[i][j];
                }
                if(i==j-1)
                {
                    c[i]=1.0-(1.0-eps)*lambda[i][j];
                }
            }

            formal_J[i] = jTerm;
            x[i] = formal_J[i] - Jold[i];
        }
        
        TridiagonalSoln(y, a, b, c, x);
        
        for (int i=0; i < nZones; i++)
        {
            nuModel.J[i] = y[i] - Jold[i];
        }        
    }

    //SMS moronic attempts at coding
    void ALIcalcS(NuModel &nuModel)
    {
        const int &nZones = nuModel.params.nZones;
        const double &eps = nuModel.params.eps;

        const vector<vector<double>> &lambda = nuModel.lambda;
        std::vector<double> Sold = nuModel.S;
        std::vector<double> a(nZones-1,0.0);
        std::vector<double> b(nZones,0.0);
        std::vector<double> c(nZones-1,0.0);
        std::vector<double> x(nZones,0.0);
        std::vector<double> y(nZones,0.0);

        for (int i=0; i < nZones; i++)
        {
            for (int j=0; j < nZones; j++)
            {
                if(i==j)
                {
                    b[i]=1.0-(1.0-eps)*lambda[i][j];
                }
                if(i==j+1)
                {
                    a[i]=1.0-(1.0-eps)*lambda[i][j];
                }
                if(i==j-1)
                {
                    c[i]=1.0-(1.0-eps)*lambda[i][j];
                }
            }
            x[i] = eps * nuModel.B[i];
        }

        TridiagonalSoln(y, a, b, c, x);

        for (int i=0; i < nZones; i++)
        {
            nuModel.S[i] = y[i] - Sold[i];
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
        ALIcalcS(nuModel);
    }
}