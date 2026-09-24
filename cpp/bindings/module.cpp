#include <pybind11/pybind11.h>

// Build/import smoke surface only. Numerical kernels arrive with oracle tests.
PYBIND11_MODULE(_core, module) {
    module.doc() = "EXCAVATOR2 numerical extension scaffold; no algorithms implemented";
    module.attr("implementation_status") = "scaffold";
}
