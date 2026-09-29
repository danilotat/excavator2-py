#include "excavator2/fastcall.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace excavator2::fastcall {
namespace {
constexpr double log_root_2pi = 0.91893853320467274178;
constexpr double root_pi_over_two = 1.25331413731550025121;
constexpr double min_deviation = 0.001;
constexpr double infinity = std::numeric_limits<double>::infinity();
constexpr double nodes[16] = {
    -0.98940093499164994,
    -0.9445750230732326,
    -0.86563120238783164,
    -0.75540440835500311,
    -0.61787624440264377,
    -0.45801677765722737,
    -0.28160355077925892,
    -0.095012509837637441,
    0.095012509837637441,
    0.28160355077925892,
    0.45801677765722737,
    0.61787624440264377,
    0.75540440835500311,
    0.86563120238783164,
    0.9445750230732326,
    0.98940093499164994
};
constexpr double weights[16] = {
    0.027152459411754055,
    0.062253523938647609,
    0.095158511682492994,
    0.12462897125553395,
    0.14959598881657665,
    0.16915651939500262,
    0.18260341504492361,
    0.18945061045506859,
    0.18945061045506859,
    0.18260341504492361,
    0.16915651939500262,
    0.14959598881657665,
    0.12462897125553395,
    0.095158511682492994,
    0.062253523938647609,
    0.027152459411754055
};

double log_mills(double x) {
    if (x < 20) {
        return std::log(root_pi_over_two * std::erfc(x / std::sqrt(2.0))) + 0.5 * x * x;
    }
    const double inverse_square = (1 / x) * (1 / x);
    double sum = 1, term = 1;
    for (int k = 1; k < 100; ++k) {
        term *= -(2 * k - 1) * inverse_square;
        sum += term;
        if (std::abs(term) < std::numeric_limits<double>::epsilon() * std::abs(sum)) break;
    }
    return -std::log(x) + std::log(sum);
}

double log_scaled_mass(double a, double b) {
    if (!std::isfinite(a) || !std::isfinite(b) || !(a < b)) {
        throw NumericalError("unrepresentable truncated interval");
    }
    if (b <= 0) {
        const double old_a = a;
        a = -b;
        b = -old_a;
    }
    const double width = b - a;
    if (width * std::max({1.0, std::abs(a), std::abs(b)}) < 1) {
        const double anchor = std::max(0.0, a);
        double sum = 0;
        for (int k = 0; k < 16; ++k) {
            const double offset = width / 2 * (nodes[k] + 1) + (a - anchor);
            sum += weights[k] * std::exp(-anchor * offset - 0.5 * offset * offset);
        }
        return std::log(width / 2) + std::log(sum);
    }
    if (a < 0) {
        return std::log(root_pi_over_two *
                        (std::erf(b / std::sqrt(2.0)) - std::erf(a / std::sqrt(2.0))));
    }
    const double left = log_mills(a), right = log_mills(b);
    return left + std::log(-std::expm1(-0.5 * (b - a) * (b + a) + right - left));
}

double normalize_row(double* row) {
    double largest = -infinity;
    for (std::size_t j = 0; j < states; ++j) {
        if (std::isnan(row[j]) || row[j] == infinity) {
            throw NumericalError("observation has no finite probability within model support");
        }
        largest = std::max(largest, row[j]);
    }
    if (!std::isfinite(largest)) {
        throw NumericalError("observation has no finite probability within model support");
    }
    double total = 0;
    for (std::size_t j = 0; j < states; ++j) {
        row[j] = std::exp(row[j] - largest);
        total += row[j];
    }
    for (std::size_t j = 0; j < states; ++j) row[j] /= total;
    return largest + std::log(total);
}

double second_moment(double sd, double lower, double upper) {
    if (sd == infinity) return (lower * lower + lower * upper + upper * upper) / 3;
    const double a = lower / sd, b = upper / sd;
    if (std::max(std::abs(a), std::abs(b)) <= 1 ||
        (b - a) * std::max({1.0, std::abs(a), std::abs(b)}) < 1) {
        const double mid = lower + (upper - lower) / 2, half = (upper - lower) / 2;
        double mass = 0, moment = 0;
        for (int k = 0; k < 16; ++k) {
            const double offset = half * nodes[k], z = offset / sd;
            const double weight = weights[k] * std::exp(-(mid / sd) * z - 0.5 * z * z);
            mass += weight;
            moment += weight * (mid + offset) * (mid + offset);
        }
        return moment / mass;
    }
    const double mass = log_scaled_mass(a, b), anchor = std::clamp(0.0, a, b);
    const double pa = std::exp(-0.5 * (a - anchor) * (a + anchor) - mass);
    const double pb = std::exp(-0.5 * (b - anchor) * (b + anchor) - mass);
    return sd * sd * (1 + a * pa - b * pb);
}

double fit_deviation(double moment, double lower, double upper) {
    if (moment >= second_moment(infinity, lower, upper)) return infinity;
    if (moment <= second_moment(min_deviation, lower, upper)) return min_deviation;
    double lo = std::log(min_deviation);
    double hi = std::log(std::max({1.0, std::abs(lower), std::abs(upper)}));
    int count = 0;
    while (second_moment(std::exp(hi), lower, upper) < moment) {
        if (++count == 128) throw NumericalError("could not bracket truncated deviation");
        hi += std::log(2.0);
    }
    for (int k = 0; k < 64; ++k) {
        const double mid = (lo + hi) / 2;
        if (second_moment(std::exp(mid), lower, upper) < moment) lo = mid;
        else hi = mid;
    }
    return std::exp((lo + hi) / 2);
}

double objective(double moment, double sd, double lower, double upper) {
    if (sd == infinity) return -std::log(upper - lower);
    const double anchor = std::clamp(0.0, lower, upper);
    return -std::log(sd) - log_scaled_mass(lower / sd, upper / sd)
           - 0.5 * (moment - anchor * anchor) / sd / sd;
}
}  // namespace

void posterior(std::size_t n, const double* values, const double* means,
               const double* deviations, const double* priors, double* output) {
    for (std::size_t i = 0; i < n; ++i) {
        double* row = output + states * i;
        for (std::size_t j = 0; j < states; ++j) {
            const double z = (values[i] - means[j]) / deviations[j];
            row[j] = priors[j] > 0 ? std::log(priors[j]) - std::log(deviations[j])
                     - log_root_2pi - 0.5 * z * z : -infinity;
        }
        normalize_row(row);
    }
}

double expectation(std::size_t n, const double* values, const double* means,
                   const double* deviations, const double* priors,
                   const double* bounds, double* output) {
    double base[states];
    for (std::size_t j = 0; j < states; ++j) {
        base[j] = -infinity;
        bool occupied = false;
        for (std::size_t i = 0; i < n; ++i) {
            if (values[i] >= bounds[2*j] && values[i] <= bounds[2*j+1]) {
                occupied = true;
                break;
            }
        }
        if (priors[j] <= 0 || !occupied) continue;
        if (deviations[j] == infinity) {
            base[j] = std::log(priors[j]) - std::log(bounds[2*j+1] - bounds[2*j]);
        } else {
            base[j] = std::log(priors[j]) - std::log(deviations[j])
                      - log_scaled_mass((bounds[2*j] - means[j]) / deviations[j],
                                        (bounds[2*j+1] - means[j]) / deviations[j]);
        }
    }
    double likelihood = 0;
    for (std::size_t i = 0; i < n; ++i) {
        double* row = output + states * i;
        for (std::size_t j = 0; j < states; ++j) {
            row[j] = -infinity;
            if (values[i] >= bounds[2*j] && values[i] <= bounds[2*j+1] && base[j] != -infinity) {
                const double anchor = std::clamp(means[j], bounds[2*j], bounds[2*j+1]);
                const double difference = (values[i] - anchor) / deviations[j];
                row[j] = base[j] - 0.5 * difference *
                         (difference + 2 * (anchor - means[j]) / deviations[j]);
            }
        }
        likelihood += normalize_row(row);
    }
    if (!std::isfinite(likelihood)) throw NumericalError("non-finite mixture log-likelihood");
    return likelihood;
}

void maximization(std::size_t n, const double* values, const double* responsibilities,
                  const double* means, const double* deviations, const double* bounds,
                  double* new_deviations, double* new_priors) {
    double total = 0;
    for (std::size_t j = 0; j < states; ++j) {
        if (!(deviations[j] >= min_deviation)) {
            throw std::invalid_argument("fitted deviations must be at least 0.001");
        }
        if (!(bounds[2*j] < bounds[2*j+1])) {
            throw std::invalid_argument("bounds must have positive width");
        }
        double mass = 0, squared_error = 0;
        for (std::size_t i = 0; i < n; ++i) {
            const double weight = responsibilities[states * i + j];
            if (weight < 0 || !std::isfinite(weight)) {
                throw std::invalid_argument("invalid responsibilities");
            }
            if (weight == 0) continue;
            if (values[i] < bounds[2*j] || values[i] > bounds[2*j+1]) {
                throw std::invalid_argument("responsibilities outside component support");
            }
            const double difference = values[i] - means[j];
            mass += weight;
            squared_error += weight * difference * difference;
        }
        new_deviations[j] = deviations[j];
        new_priors[j] = mass;
        total += mass;
        if (mass == 0) continue;
        const double lower = bounds[2*j] - means[j], upper = bounds[2*j+1] - means[j];
        const double moment = squared_error / mass;
        const double candidate = fit_deviation(moment, lower, upper);
        const double old_q = objective(moment, deviations[j], lower, upper);
        const double new_q = objective(moment, candidate, lower, upper);
        if (!std::isfinite(old_q) || !std::isfinite(new_q)) {
            throw NumericalError("non-finite truncated objective");
        }
        if (new_q >= old_q) new_deviations[j] = candidate;
    }
    if (!(total > 0) || !std::isfinite(total)) {
        throw std::invalid_argument("responsibilities must have finite positive total mass");
    }
    for (std::size_t j = 0; j < states; ++j) new_priors[j] /= total;
}
}  // namespace excavator2::fastcall
