#pragma once

#include <cstddef>

namespace excavator2::fastcall {
constexpr std::size_t states = 5;

// All arrays are contiguous binary64. Matrices are row-major: (n, 5),
// except bounds (5, 2). Inputs and outputs must not alias; n must be positive.
// The Python binding validates dimensions and parameter domains.
void posterior(std::size_t n, const double* values, const double* means,
               const double* deviations, const double* priors, double* output);
void expectation(std::size_t n, const double* values, const double* means,
                 const double* deviations, const double* priors,
                 const double* bounds, double* output);
void maximization(std::size_t n, const double* values, const double* responsibilities,
                  const double* means, const double* deviations,
                  double* new_deviations, double* new_priors);
}  // namespace excavator2::fastcall
