#include <iostream>
#include <vector>

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
                         const vector<double> &a, const vector<double> &b,
                         const vector<double> &c, const vector<double> &x,
                         vector<double> &cHelper, vector<double> &dHelper)
    {
        const int N = x.size();

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
}