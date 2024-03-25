#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "EnhancementFactor.hpp"
#include "System.hpp"

namespace py = pybind11;
using namespace LED;

PYBIND11_MODULE(LED, m) {
  py::class_<System>(m, "System")
      .def(py::init([](const std::string &filepath) {
             return new System(filepath);
           }),
           py::arg("filepath") = std::string(""),
           py::return_value_policy::take_ownership)
      .def_readonly("name", &System::name);

  py::class_<EnhancementFactor>(m, "EnhancementFactor")
      .def(py::init([](int factorization_order, bool use_tau, bool use_nrg) {
             return new EnhancementFactor(factorization_order, use_tau,
                                          use_nrg);
           }),
           py::arg("factorization_order") = 1, py::arg("use_tau") = false,
           py::arg("use_nrg") = false, py::return_value_policy::take_ownership)
      .def_property("coefficients", &EnhancementFactor::get_coefficients,
                    &EnhancementFactor::set_coefficients)
      .def_property_readonly("N_free_params",
                             &EnhancementFactor::get_N_free_parameters)
      .def("compute_error",
           py::overload_cast<const System &, double>(
               &EnhancementFactor::compute_error),
           py::arg("system"), py::arg("omega") = 0.5)
      .def("compute_error",
           py::overload_cast<const std::vector<System> &, double>(
               &EnhancementFactor::compute_error),
           py::arg("systems_list"), py::arg("omega") = 0.5)
      .def("compute_total_energy", &EnhancementFactor::compute_total_energy);
}
