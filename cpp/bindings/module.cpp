#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <cmath>
#include <cstdint>
#include <string>

#include "excavator2/fastcall.hpp"

namespace py = pybind11;
namespace fc = excavator2::fastcall;
// noconvert() below rejects wrong dtypes/strides instead of hiding copies.
using Array = py::array_t<double, py::array::c_style>;

namespace {
void check_array(const Array& array, const char* name, py::ssize_t rows,
                 py::ssize_t columns = 0) {
    const bool shape_ok = columns == 0
        ? array.ndim() == 1 && array.shape(0) == rows
        : array.ndim() == 2 && array.shape(0) == rows && array.shape(1) == columns;
    if (!shape_ok) throw py::value_error(std::string(name) + " has an invalid shape");
    if (reinterpret_cast<std::uintptr_t>(array.data()) % alignof(double) != 0) {
        throw py::value_error(std::string(name) + " must be aligned for binary64");
    }
    for (py::ssize_t i = 0; i < array.size(); ++i) {
        if (!std::isfinite(array.data()[i])) {
            throw py::value_error(std::string(name) + " must contain finite values");
        }
    }
}

py::ssize_t check_inputs(const Array& values, const Array& means, const Array& deviations) {
    if (values.ndim() != 1 || values.size() == 0) {
        throw py::value_error("values must be a nonempty vector");
    }
    const auto n = values.size();
    check_array(values, "values", n);
    check_array(means, "means", fc::states);
    check_array(deviations, "deviations", fc::states);
    for (py::ssize_t j = 0; j < static_cast<py::ssize_t>(fc::states); ++j) {
        if (deviations.data()[j] <= 0) throw py::value_error("deviations must be positive");
    }
    return n;
}

void check_nonnegative(const Array& array, const char* name) {
    for (py::ssize_t i = 0; i < array.size(); ++i) {
        if (array.data()[i] < 0) throw py::value_error(std::string(name) + " must be nonnegative");
    }
}

Array posterior(const Array& values, const Array& means, const Array& deviations,
                const Array& priors) {
    const auto n = check_inputs(values, means, deviations);
    check_array(priors, "priors", fc::states);
    check_nonnegative(priors, "priors");
    Array result({n, static_cast<py::ssize_t>(fc::states)});
    const auto* x = values.data(); const auto* mu = means.data();
    const auto* sd = deviations.data(); const auto* p = priors.data();
    auto* output = result.mutable_data();
    {
        py::gil_scoped_release release;
        fc::posterior(n, x, mu, sd, p, output);
    }
    return result;
}

Array expectation(const Array& values, const Array& means, const Array& deviations,
                  const Array& priors, const Array& bounds) {
    const auto n = check_inputs(values, means, deviations);
    check_array(priors, "priors", fc::states);
    check_nonnegative(priors, "priors");
    check_array(bounds, "bounds", fc::states, 2);
    for (std::size_t j = 0; j < fc::states; ++j) {
        if (bounds.data()[2*j] > bounds.data()[2*j+1]) throw py::value_error("unordered bounds");
    }
    Array result({n, static_cast<py::ssize_t>(fc::states)});
    const auto* x = values.data(); const auto* mu = means.data();
    const auto* sd = deviations.data(); const auto* p = priors.data();
    const auto* limits = bounds.data(); auto* output = result.mutable_data();
    {
        py::gil_scoped_release release;
        fc::expectation(n, x, mu, sd, p, limits, output);
    }
    return result;
}

py::tuple maximization(const Array& values, const Array& responsibilities,
                       const Array& means, const Array& deviations) {
    const auto n = check_inputs(values, means, deviations);
    check_array(responsibilities, "responsibilities", n, fc::states);
    check_nonnegative(responsibilities, "responsibilities");
    Array sd_out(fc::states), priors_out(fc::states);
    const auto* x = values.data(); const auto* weights = responsibilities.data();
    const auto* mu = means.data(); const auto* sd = deviations.data();
    auto* new_sd = sd_out.mutable_data(); auto* new_priors = priors_out.mutable_data();
    {
        py::gil_scoped_release release;
        fc::maximization(n, x, weights, mu, sd, new_sd, new_priors);
    }
    return py::make_tuple(sd_out, priors_out);
}
}  // namespace

PYBIND11_MODULE(_core, module) {
    module.doc() = "Scalar FastCall numerical kernels; algorithm control remains in Python";
    module.def("fastcall_posterior", &posterior, py::arg("values").noconvert(),
               py::arg("means").noconvert(), py::arg("deviations").noconvert(),
               py::arg("priors").noconvert());
    module.def("fastcall_expectation", &expectation, py::arg("values").noconvert(),
               py::arg("means").noconvert(), py::arg("deviations").noconvert(),
               py::arg("priors").noconvert(), py::arg("bounds").noconvert());
    module.def("fastcall_maximization", &maximization, py::arg("values").noconvert(),
               py::arg("responsibilities").noconvert(), py::arg("means").noconvert(),
               py::arg("deviations").noconvert());
}
