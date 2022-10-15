#include <iostream>
#include <vector>

using namespace std;

namespace myLib
{
    const double trapezoidal(const vector<double> X,
                             const vector<double> integrand)
    {
        const int nPoints = X.size();
        double integral = 0.0;

        for (int i = 1; i < nPoints; i++)
        {
            integral += (integrand[i - 1] + integrand[i]) * (X[i] - X[i - 1]);
        }

        integral *= 0.5;
        return integral;
    }
}