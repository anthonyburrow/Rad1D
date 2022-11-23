#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "RadModel.hpp"
#include "../tests/test_cpp/test_main.hpp"

using namespace std;
namespace py = pybind11;

PYBIND11_MODULE(Rad1D, module_handle) {
    module_handle.doc() = "Radiative transfer model.";
    // RadModel Class
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
    // Testing
    module_handle.def("test_all", &myLib::testAll);
}
