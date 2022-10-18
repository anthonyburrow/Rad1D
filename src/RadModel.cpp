#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "RadModel.hpp"
#include "NuModel.hpp"
#include "io.hpp"
#include "initialize.hpp"
#include "constants.hpp"
#include "lambda.hpp"

using namespace std;
namespace py = pybind11;

namespace myLib
{
    RadModel::RadModel(const py::dict &dictParams)
      : params(dictToParams(dictParams)),
        // Allocate vectors
        quadMu(params.nQuad, zero),
        quadW(params.nQuad, zero),
        spectrum(params.nWave, vector<double>(2, zero)),
        tau(params.nZones, zero),
        T(params.nZones, zero),
        lambda(params.nZones, vector<double>(params.nZones, zero))
        // lambdaA(params.nZones, zero),
        // lambdaB(params.nZones, zero),
        // lambdaC(params.nZones, zero)
    {
        initialize(*this);
        initLambda();
    }

    void RadModel::initLambda()
    {
        cout << "Calculating Lambda matrix..." << endl;
        auto t0 = chrono::high_resolution_clock::now();
        chrono::duration<double> sec;

        // calcLambda(*this, lambdaA, lambdaB, lambdaC);
        calcLambda(*this, lambda);

        auto t1 = chrono::high_resolution_clock::now();
        sec = t1 - t0;
        // cout << "Calculated Lambda matrix in "<< fixed << setprecision(3)
        //      << sec.count() << " sec" << endl;
    }

    vector<vector<double>> RadModel::genSpectrum(bool normalize)
    {
        cout << "Generating spectrum..." << endl;
        auto t0 = chrono::high_resolution_clock::now();
        chrono::duration<double> sec;

        // Calc flux at each wavepoint - parallelize this in the future?
        for (int i=0; i < params.nWave; i++)
        {
            NuModel nuModel = NuModel(spectrum[i][0], *this);
            spectrum[i][1] = nuModel.calcFlux();
        }

        if (normalize) { normalizeFlux(); }
        else { scaleFlux(); }

        auto t1 = chrono::high_resolution_clock::now();
        sec = t1 - t0;
        cout << "Generated spectrum in "<< fixed << setprecision(3)
             << sec.count() << " sec" << endl;

        return spectrum;
    }

    void RadModel::normalizeFlux()
    {
        double maxFlux = (*max_element(begin(spectrum), end(spectrum),
                          [](auto &a, auto &b){ return a[1] < b[1]; }))[1];

        for (int i=0; i < params.nWave; i++)
        {
            spectrum[i][1] /= maxFlux;
        }
    }

    void RadModel::scaleFlux()
    {
        // Rescale with constants from Planck function
        const double scale = 2.0 * hc * c * 1e8;

        for (int i=0; i < params.nWave; i++)
        {
            spectrum[i][1] *= scale;
        }
    }

    vector<vector<double>> RadModel::convergenceTest(const double &lam)
    {
        vector<vector<double>> results(params.maxIter,
                                       vector<double>(params.nZones));

        NuModel nuModel = NuModel(lam, *this);

        for (int i=0; i < params.maxIter; i++)
        {
            nuModel.iterate();

            for (int j=0; j < params.nZones; j++)
            {
                results[i][j] = nuModel.J[j] / nuModel.B[j];
            }
        }

        return results;
    }

    // Properties
    std::vector<double> RadModel::getTau() { return tau; }
    std::vector<double> RadModel::getT() { return T; }
}

PYBIND11_MODULE(Rad1D, module_handle) {
    module_handle.doc() = "Radiative transfer model.";
    py::class_<myLib::RadModel>(module_handle, "RadModel")
        .def(py::init<const py::dict &>()
        )
        .def("gen_spectrum", [](myLib::RadModel &self, bool normalize) {
            py::array out = py::cast(self.genSpectrum(normalize));
            return out;
        })
        .def("convergence_test", [](myLib::RadModel &self, double &lam) {
            py::array out = py::cast(self.convergenceTest(lam));
            return out;
        })
        // Properties
        .def_property_readonly("tau", [](myLib::RadModel &self) {
            py::array out = py::cast(self.getTau()); return out;
        })
        .def_property_readonly("T", [](myLib::RadModel &self) {
            py::array out = py::cast(self.getT()); return out;
        })
    ;
}
