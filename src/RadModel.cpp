#include <iostream>
#include <string>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "io.hpp"
#include "RadModel.hpp"
#include "NuModel.hpp"
#include "initialize.hpp"
#include "constants.hpp"
#include "util.hpp"
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
        tauCont(params.nZones, zero),
        T(params.nZones, zero)
    {
        readLineList(*this);
        initialize(*this);
    }

    vector<vector<double>> RadModel::genSpectrum(bool normalize)
    {
        ostringstream output1;
        output1 << "Generating spectrum...";
        log(output1);

        auto t0 = chrono::high_resolution_clock::now();
        chrono::duration<double> sec;

        // Populate wavelengths
        vector<vector<double>> spectrum = initSpectrum(*this);
        const int nWave = static_cast<int>(spectrum.size());

        // Calc flux at each wavepoint - parallelize this in the future?
        for (int i=0; i < nWave; i++)
        {
            NuModel nuModel = NuModel(spectrum[i][0], *this);
            spectrum[i][1] = nuModel.calcFlux();
        }

        if (normalize) { normalizeSpec(spectrum); }

        auto t1 = chrono::high_resolution_clock::now();
        sec = t1 - t0;

        ostringstream output2;
        output2 << "  Finished in "<< fixed << setprecision(3)
               << sec.count() << " sec";
        log(output2);

        return spectrum;
    }

    vector<vector<double>> RadModel::convergenceTest(const double &lam)
    {
        ostringstream output1;
        output1 << "Performing convergence test at " << lam << " A...";
        log(output1);

        auto t0 = chrono::high_resolution_clock::now();
        chrono::duration<double> sec;

        vector<vector<double>> results(params.maxIter,
                                       vector<double>(params.nZones));

        NuModel nuModel = NuModel(lam, *this);

        for (int j=0; j < params.nZones; j++)
        {
            results[0][j] = nuModel.S[j] / nuModel.B[j];
        }

        for (int i=1; i < params.maxIter; i++)
        {
            nuModel.iterate();

            for (int j=0; j < params.nZones; j++)
            {
                results[i][j] = nuModel.S[j] / nuModel.B[j];
            }
        }

        auto t1 = chrono::high_resolution_clock::now();
        sec = t1 - t0;

        ostringstream output2;
        output2 << "  Finished in "<< fixed << setprecision(3)
               << sec.count() << " sec";
        log(output2);

        return results;
    }

    void RadModel::log(const ostringstream &output) const
    {
        if (!params.verbose) { return; }

        // Create log file and send to file too?
        cout << output.str() << endl;
    }

    // Properties
    std::vector<double> RadModel::getTau() { return tauCont; }
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
        },
            py::arg("normalize") = false)
        .def("convergence_test", [](myLib::RadModel &self, double &lam) {
            py::array out = py::cast(self.convergenceTest(lam));
            return out;
        },
            py::arg("lam") = 5000.0)
        // Properties
        .def_property_readonly("tau", [](myLib::RadModel &self) {
            py::array out = py::cast(self.getTau()); return out;
        })
        .def_property_readonly("T", [](myLib::RadModel &self) {
            py::array out = py::cast(self.getT()); return out;
        })
    ;
}
