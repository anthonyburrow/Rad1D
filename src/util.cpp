#include <iostream>
#include <vector>
#include <limits>
#include <cmath>

// Not sure why I need this here...
#include <pybind11/numpy.h>

#include "util.hpp"

using namespace std;

namespace myLib
{
    const double trapezoidal(const vector<double> X,
                             const vector<double> integrand)
    {
        double integral = 0.0;

        for (int i = 1; i < X.size(); i++)
        {
            integral += (integrand[i - 1] + integrand[i]) * (X[i] - X[i - 1]);
        }

        integral *= 0.5;
        return integral;
    }

    void normalizeSpec(vector<vector<double>> &spectrum)
    {
        double maxFlux = (*max_element(begin(spectrum), end(spectrum),
                          [](auto &a, auto &b){ return a[1] < b[1]; }))[1];

        for (vector<double> &point : spectrum)
        {
            point[1] /= maxFlux;
        }
    }

    const int closestIndex(const vector<double> &vec, const double value) {
        auto const it = lower_bound(vec.begin(), vec.end(), value);
        const int ind = it - vec.begin();

        return ind;
    }

    void TridiagonalSoln(vector<double> &y,
                         const double &eps, 
                         const std::vector<std::vector<double>> &lambda, 
                         const vector<double> &x
                         )
    {
        
        std::vector<double> cHelper;
        std::vector<double> dHelper;
        
        const int N = x.size();
        std::vector<double> a(N-1,0.0);
        std::vector<double> b(N,0.0);
        std::vector<double> c(N-1,0.0);

        // a is lower diagonal, c is upper diagonal
        for (int i=0; i < N; i++)
        {
            b[i]=1.0-(1.0-eps)*lambda[i][i];
        }
        // this should be indexing the correct index now
        for(int i = 0; i < N-1; i++)
        {
            a[i]= -(1.0-eps)*lambda[i+1][i];
        }
        for(int i = 0; i < N-1; i++)
        {
            c[i]= -(1.0-eps)*lambda[i][i+1];
        }

        cHelper.resize(N);
        dHelper.resize(N);

        cHelper[0] = c[0] / b[0];
        dHelper[0] = x[0] / b[0];

        // Forward sweep
        double m;
        for (int i = 1; i < N; i++) {
            m = 1.0 / (b[i] - a[i] * cHelper[i - 1]);
            cHelper[i] = c[i] * m;
            dHelper[i] = (x[i] - a[i] * dHelper[i - 1]) * m;
        }

        // Reverse sweep
        for (int i = N - 1; i-- > 0; ) {
            y[i] = dHelper[i] - cHelper[i] * x[i + 1];
        }
        
    }

    const vector<double> expn(const int &n, const double &x)
    {
        vector<double> En(n + 1, 0.0);

        // Re-write of specfun.f ENXB subroutine
        if (x == 0.0)
        {
            En[0] = numeric_limits<double>::quiet_NaN();
            En[1] = numeric_limits<double>::quiet_NaN();

            for (int k = 2; k < n + 1; k++)
            {
                En[k] = 1.0 / (k - 1.0);
            }

            return En;
        }

        En[0] = exp(-x) / x;

        if (x <= 1.0)
        {
            double S0, S, RP, PS, ENS, R;

            S0 = 0.0;
            for (int l = 1; l < n + 1; l++)
            {
                RP = 1.0;
                for (int j = 1; j < l; j++) { RP = -RP * x / j; }
                PS = -0.5772156649015328;
                for (int m = 1; m < l; m++) { PS += 1.0 / m; }
                ENS = RP * (-log(x) + PS);
                S = 0.0;
                for (int m = 0; m < 21; m++)
                {
                    if (m == l - 1) { continue; }
                    R = 1.0;
                    for (int j = 1; j < m + 1; j++) { R *= -x / j; }
                    S += R / (m - l + 1.0);
                    if (abs(S - S0) < abs(S) * 1e-15) { break; }
                    S0 = S;
                }
                En[l] = ENS - S;
            }

            return En;
        }

        double M, T0, T;

        M = 15 + int(100.0 / x);
        for (int l = 1; l < n + 1; l++)
        {
            T0 = 0.0;
            for (int k = M; k > 0; k--)
            {
                T0 = (l + k - 1.0) / (1.0 + k / (x + T0));
            }
            T = 1.0 / (x + T0);
            En[l] = exp(-x) * T;
        }

        return En;
    }
}