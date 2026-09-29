#pragma once
#include <cstddef>
#include <cstdint>

namespace excavator2::hslm {
// Single-profile binary64 buffers; n,k > 0. eta has n-1 entries.
// Outputs: one-based path[n], optional transitions[n-1,k,k], emissions[n,k],
// scores[n,k], predecessors[n,k]. All matrices use row-major Python orientation.
// Optional trace pointers may be null. Inputs/outputs must not alias.
void segment(std::size_t n, std::size_t k, const double* values, const double* means,
             double mi, double smu, double sepsilon, const double* eta,
             const double* initial, std::int32_t* path, double* transitions,
             double* emissions, double* scores, std::int32_t* predecessors);
}
