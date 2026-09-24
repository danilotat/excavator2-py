#include "excavator2/fastcall.hpp"

#include <cmath>
#include <limits>

namespace excavator2::fastcall {
namespace {
constexpr double pi = 3.141592653589793238462643383279502884;

// Ordinary probability space is deliberate: preserve legacy underflow/fallback.
double density(double value, double mean, double deviation) {
    const double z = (value - mean) / deviation;
    return std::exp(-0.5 * (z * z)) / std::sqrt(2.0 * pi) / deviation;
}

double distribution(double z) {
    return 0.5 * std::erfc(-z / std::sqrt(2.0));
}

void normalize_row(double* row, double value, const double* means) {
    double total = 0.0;
    for (std::size_t j = 0; j < states; ++j) {
        // R rowSums(na.rm=TRUE) ignores NaNs without deleting their numerators.
        if (!std::isnan(row[j])) total += row[j];
    }
    if (total == 0.0) {
        std::size_t nearest = 0;
        double distance = std::abs(means[0] - value);
        for (std::size_t j = 1; j < states; ++j) {
            const double candidate = std::abs(means[j] - value);
            if (candidate < distance) {
                nearest = j;
                distance = candidate;
            }
        }
        row[nearest] = 1.0;
        total = 0.0;
        for (std::size_t j = 0; j < states; ++j) {
            if (!std::isnan(row[j])) total += row[j];
        }
    }
    for (std::size_t j = 0; j < states; ++j) row[j] /= total;
}
}  // namespace

void posterior(std::size_t n, const double* values, const double* means,
               const double* deviations, const double* priors, double* output) {
    for (std::size_t i = 0; i < n; ++i) {
        double* row = output + states * i;
        for (std::size_t j = 0; j < states; ++j) {
            row[j] = priors[j] * density(values[i], means[j], deviations[j]);
        }
        normalize_row(row, values[i], means);
    }
}

void expectation(std::size_t n, const double* values, const double* means,
                 const double* deviations, const double* priors,
                 const double* bounds, double* output) {
    double denominator[states] = {};
    bool occupied[states] = {};
    for (std::size_t j = 0; j < states; ++j) {
        const double lower = bounds[2 * j];
        const double upper = bounds[2 * j + 1];
        for (std::size_t i = 0; i < n; ++i) {
            if (values[i] >= lower && values[i] <= upper) {
                occupied[j] = true;
                break;
            }
        }
        if (occupied[j]) {
            denominator[j] = distribution((upper - means[j]) / deviations[j])
                           - distribution((lower - means[j]) / deviations[j]);
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        double* row = output + states * i;
        for (std::size_t j = 0; j < states; ++j) {
            double term = 0.0;
            if (occupied[j]) {
                const bool inside = values[i] >= bounds[2 * j] && values[i] <= bounds[2 * j + 1];
                // Do not short-circuit outside values: 0/0 must remain NaN.
                term = density(values[i], means[j], deviations[j]) * inside / denominator[j];
                if (term == std::numeric_limits<double>::infinity()) term = 100.0;
            }
            row[j] = term * priors[j];
        }
        normalize_row(row, values[i], means);
    }
}

void maximization(std::size_t n, const double* values, const double* responsibilities,
                  const double* means, const double* deviations,
                  double* new_deviations, double* new_priors) {
    double prior_total = 0.0;
    for (std::size_t j = 0; j < states; ++j) {
        double mass = 0.0;
        double squared_error = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double weight = responsibilities[states * i + j];
            const double difference = values[i] - means[j];
            mass += weight;
            squared_error += weight * (difference * difference);
        }
        new_deviations[j] = deviations[j];
        if (mass != 0.0) {
            const double candidate = std::sqrt(squared_error / mass);
            if (!(candidate < 1e-100)) new_deviations[j] = candidate;
            new_priors[j] = mass / static_cast<double>(n);
        } else {
            new_priors[j] = 1e-6;
        }
        prior_total += new_priors[j];
    }
    for (std::size_t j = 0; j < states; ++j) new_priors[j] /= prior_total;
}
}  // namespace excavator2::fastcall
