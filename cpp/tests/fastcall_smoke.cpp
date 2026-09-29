// Standalone memory-safety smoke test; numerical oracle tests live in Python.
#include "excavator2/fastcall.hpp"
#include <array>
#include <cassert>
#include <cmath>
#include <vector>

int main() {
    using namespace excavator2::fastcall;
    const std::array<double, 5> means{-3, -1, 0, 0.58, 1};
    const std::array<double, 5> sd{0.1, 0.1, 0.1, 0.1, 0.1};
    const std::array<double, 5> priors{0.05, 0.1, 0.7, 0.1, 0.05};
    const std::array<double, 10> bounds{-20, -1.3, -1.3, -0.5, -0.5, 0.35, 0.35, 0.9, 0.9, 20};
    for (std::size_t n : {1, 17, 10003}) {
        std::vector<double> values(n), weights(n * states), probabilities(n * states);
        for (std::size_t i = 0; i < n; ++i) values[i] = means[i % states];
        expectation(n, values.data(), means.data(), sd.data(), priors.data(), bounds.data(), weights.data());
        std::array<double, 5> next_sd{}, next_priors{};
        maximization(n, values.data(), weights.data(), means.data(), sd.data(), next_sd.data(), next_priors.data());
        posterior(n, values.data(), means.data(), next_sd.data(), next_priors.data(), probabilities.data());
        for (std::size_t i = 0; i < n; ++i) {
            double sum = 0;
            for (std::size_t j = 0; j < states; ++j) sum += probabilities[states * i + j];
            assert(std::isfinite(sum) && std::abs(sum - 1) < 1e-14);
        }
    }
}
