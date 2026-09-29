#pragma once

#include <cstddef>
#include <stdexcept>

namespace excavator2::fastcall {
constexpr std::size_t states = 5;
class NumericalError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// All arrays are contiguous binary64. Matrices are row-major: (n, 5),
// except bounds (5, 2). Inputs and outputs must not alias; n must be positive.
// The Python binding validates dimensions and parameter domains.
void posterior(std::size_t n, const double* values, const double* means,
               const double* deviations, const double* priors, double* output);
double expectation(std::size_t n, const double* values, const double* means,
                 const double* deviations, const double* priors,
                 const double* bounds, double* output);
void maximization(std::size_t n, const double* values, const double* responsibilities,
                  const double* means, const double* deviations, const double* bounds,
                  double* new_deviations, double* new_priors);
}  // namespace excavator2::fastcall
