#include "excavator/fastcall.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <vector>

namespace excavator {
namespace fastcall {

// -----------------------------------------------------------------------------
// Helper math functions
// -----------------------------------------------------------------------------

static const double PI = 3.14159265358979323846;

// Standard normal CDF using error function
static double normal_cdf(double x, double mean, double sd) {
    return 0.5 * (1.0 + std::erf((x - mean) / (sd * std::sqrt(2.0))));
}

// Standard normal PDF
static double normal_pdf(double x, double mean, double sd) {
    double z = (x - mean) / sd;
    return std::exp(-0.5 * z * z) / (sd * std::sqrt(2.0 * PI));
}

// -----------------------------------------------------------------------------
// FastCall Implementation
// -----------------------------------------------------------------------------

FastCall::FastCall(const FastCallParameters& params)
    : params_(params)
{
    // Initialize with defaults
    means_ = DEFAULT_MEANS;
    priors_ = DEFAULT_PRIORS;
    sds_.fill(0.01);
}

FastCallResult FastCall::call(
    const std::vector<double>& segment_means,
    const std::vector<double>& segment_sds
) {
    FastCallResult result;
    result.success = false;

    // Input validation
    if (segment_means.empty()) {
        result.error_message = "Empty segment data";
        return result;
    }

    // 1. Initialize Starting Conditions (StartCond in R)
    initialize_start_conditions(segment_means);

    // 2. Run EM Algorithm (EMFastCall in R)
    auto [iterations, converged] = run_em(segment_means);
    
    result.iterations = iterations;
    result.converged = converged;

    // 3. Compute final posteriors
    auto posteriors = compute_posteriors(segment_means);

    // 4. Assign labels
    result.calls = assign_labels(posteriors, segment_means);

    // 5. Store fitted parameters
    result.state_means.assign(means_.begin(), means_.end());
    result.state_sds.assign(sds_.begin(), sds_.end());
    result.state_priors.assign(priors_.begin(), priors_.end());

    result.success = true;
    return result;
}

void FastCall::initialize_start_conditions(const std::vector<double>& data) {
    // Ported from StartCond in LibraryFastCall.R
    
    double thrd = params_.thrd;
    double thru = params_.thru;

    // Means are fixed
    means_ = DEFAULT_MEANS;

    // R Code 'StartCond' boundaries:
    // lvec <- c(-50, -1.5, -thrd, thru, 0.9)
    // uvec <- c(-1.5, -thrd, thru, 0.9, 50)
    lower_bounds_ = {-50.0, -1.5, -thrd, thru, 0.9};
    upper_bounds_ = {-1.5, -thrd, thru, 0.9, 50.0};

    // R Code: sdvec<-c(0.01,0.01,0.01,0.01,0.01)
    sds_.fill(0.01);

    for (int i = 0; i < N_STATES; ++i) {
        double lower = lower_bounds_[i];
        double upper = upper_bounds_[i];

        // R Code: ind<-which(mdata<=u & mdata>=l)
        std::vector<double> region_data;
        region_data.reserve(data.size());
        for (double x : data) {
            if (x >= lower && x <= upper) {
                region_data.push_back(x);
            }
        }

        // R Code: if (length(ind)>1) { sdvec[i]<-sd(mdata[ind]) }
        // R Code: if length is 0 or 1, it keeps the previous value (0.01)
        if (region_data.size() > 1) {
            double sum = std::accumulate(region_data.begin(), region_data.end(), 0.0);
            double mean = sum / region_data.size();
            double sq_sum = 0.0;
            for (double x : region_data) {
                sq_sum += (x - mean) * (x - mean);
            }
            double sd = std::sqrt(sq_sum / (region_data.size() - 1));
            
            sds_[i] = sd; // Assign calculated SD regardless of value here
        }
    }

    // R Code: sdvec[which(sdvec<0.001)]<-0.001
    // This happens AFTER all assignments
    for (int i = 0; i < N_STATES; ++i) {
        if (sds_[i] < 0.001) {
            sds_[i] = 0.001;
        }
    }
}

std::pair<int, bool> FastCall::run_em(const std::vector<double>& data) {
    // Ported from EMFastCall in LibraryFastCall.R

    // 1. StartCond is already called by the public method `call`
    
    // 2. R Code: prior<-c(0.05,0.1,0.7,0.1,0.05)
    priors_ = DEFAULT_PRIORS;

    // 3. R Code (EM Loop Boundaries re-definition):
    // lvec<-c(-20,-1.3,-thrd,thru,0.9)
    // uvec<-c(-1.3,-thrd,thru,0.9,20)
    // NOTE: This differs from Initialization boundaries (-50 vs -20, -1.5 vs -1.3)
    double thrd = params_.thrd;
    double thru = params_.thru;
    lower_bounds_ = {-20.0, -1.3, -thrd, thru, 0.9};
    upper_bounds_ = {-1.3, -thrd, thru, 0.9, 20.0};

    // 4. Initial "Likelihood" check
    // R Code: LikeliNew<-sum(PosteriorP(mdata,muvec,sdvec,prior)*prior)
    auto posteriors = compute_posteriors(data);
    double likelihood_old = compute_log_likelihood(data, posteriors);

    double threshold = params_.convergence;
    int max_iter = params_.max_iterations;

    int iter;
    bool converged = false;

    // R Code: for (i in 1:1000)
    for (iter = 0; iter < max_iter; ++iter) {
        // E-step
        auto tau = e_step(data);

        // M-step
        m_step(data, tau);

        // Compute new "Likelihood"
        posteriors = compute_posteriors(data);
        double likelihood_new = compute_log_likelihood(data, posteriors);

        // Check convergence
        // R Code: if (abs(LikeliNew-LikeliOld)<threshold)
        if (std::abs(likelihood_new - likelihood_old) < threshold) {
            converged = true;
            break;
        }

        likelihood_old = likelihood_new;
    }

    return {iter + 1, converged};
}

std::vector<std::vector<double>> FastCall::e_step(const std::vector<double>& data) {
    // Ported from EStep in LibraryFastCall.R
    
    size_t n = data.size();
    std::vector<std::vector<double>> tau(n, std::vector<double>(N_STATES, 0.0));

    for (size_t i = 0; i < n; ++i) {
        double x = data[i];
        std::vector<double> probs(N_STATES);
        double sum_prob = 0.0;

        for (int j = 0; j < N_STATES; ++j) {
            double lower = lower_bounds_[j];
            double upper = upper_bounds_[j];
            
            // R Code: if (sum((mdata<=u)*(mdata>=l))!=0) check
            // R gfct function uses (x<=u)*(x>=l)
            bool in_range = (x >= lower && x <= upper);

            if (in_range) {
                double pdf = truncated_gaussian_pdf(x, means_[j], sds_[j], lower, upper);
                
                // R Code: normaldataVec[which(normaldataVec==Inf)]<-100
                if (std::isinf(pdf)) {
                    pdf = 100.0;
                }
                probs[j] = priors_[j] * pdf;
            } else {
                probs[j] = 0.0;
            }
            sum_prob += probs[j];
        }

        // R Code: ind0<-which(deno==0) ... indmin<-which.min(...) ... tauMat[ind0[k],indmin]<-1
        if (sum_prob == 0.0) {
            int closest = 0;
            double min_dist = std::abs(means_[0] - x);
            for (int j = 1; j < N_STATES; ++j) {
                double dist = std::abs(means_[j] - x);
                if (dist < min_dist) {
                    min_dist = dist;
                    closest = j;
                }
            }
            tau[i][closest] = 1.0;
        } else {
            for (int j = 0; j < N_STATES; ++j) {
                tau[i][j] = probs[j] / sum_prob;
            }
        }
    }

    return tau;
}

void FastCall::m_step(
    const std::vector<double>& data,
    const std::vector<std::vector<double>>& tau
) {
    // Ported from MStep in LibraryFastCall.R
    size_t n = data.size();

    // R Code: ptmp <- colSums(taux)
    std::array<double, N_STATES> col_sums;
    col_sums.fill(0.0);
    for (size_t i = 0; i < n; ++i) {
        for (int j = 0; j < N_STATES; ++j) {
            col_sums[j] += tau[i][j];
        }
    }

    for (int j = 0; j < N_STATES; ++j) {
        if (col_sums[j] > 0) {
            // R Code: moy[j] <- muvec[j] (Means are fixed)
            
            // R Code: sqrt((taux[,j]%*%((mdata-moy[j])^2))/ptmp[j])
            double weighted_sq_sum = 0.0;
            for (size_t i = 0; i < n; ++i) {
                double diff = data[i] - means_[j];
                weighted_sq_sum += tau[i][j] * diff * diff;
            }
            double new_sd = std::sqrt(weighted_sq_sum / col_sums[j]);

            // R Code: if (... < 1e-100) { sdev[j]<-sdvec[j] } else { sdev[j]<- ... }
            if (new_sd >= 1e-100) {
                sds_[j] = new_sd;
            }

            // R Code: pnew[j]<-ptmp[j]/length(mdata)
            priors_[j] = col_sums[j] / n;
        } else {
            // R Code: else { ... pnew[j]<-1e-06 }
            priors_[j] = 1e-6;
        }
    }

    // R Code: pnew<-pnew/sum(pnew)
    double prior_sum = std::accumulate(priors_.begin(), priors_.end(), 0.0);
    for (int j = 0; j < N_STATES; ++j) {
        priors_[j] /= prior_sum;
    }
}

std::vector<std::vector<double>> FastCall::compute_posteriors(
    const std::vector<double>& data
) {
    // Ported from PosteriorP in LibraryFastCall.R
    size_t n = data.size();
    std::vector<std::vector<double>> posteriors(n, std::vector<double>(N_STATES));

    for (size_t i = 0; i < n; ++i) {
        double x = data[i];
        double sum_prob = 0.0;
        std::vector<double> probs(N_STATES);

        for (int j = 0; j < N_STATES; ++j) {
            double pdf = normal_pdf(x, means_[j], sds_[j]);
            probs[j] = priors_[j] * pdf;
            sum_prob += probs[j];
        }

        // R Code: deno <- rowSums(tauMat, na.rm = T); ind0 <- which(deno == 0) ...
        if (sum_prob == 0.0) {
            int closest = 0;
            double min_dist = std::abs(means_[0] - x);
            for (int j = 1; j < N_STATES; ++j) {
                double dist = std::abs(means_[j] - x);
                if (dist < min_dist) {
                    min_dist = dist;
                    closest = j;
                }
            }
            posteriors[i][closest] = 1.0;
        } else {
            for (int j = 0; j < N_STATES; ++j) {
                posteriors[i][j] = probs[j] / sum_prob;
            }
        }
    }
    return posteriors;
}

double FastCall::truncated_gaussian_pdf(
    double x, double mean, double sd, double lower, double upper
) {
    // Ported from gfct in LibraryFastCall.R
    // R: (dnorm(...) * (x<=u) * (x>=l)) / (pnorm(u...) - pnorm(l...))
    
    if (x < lower || x > upper) return 0.0;

    double cdf_upper = normal_cdf(upper, mean, sd);
    double cdf_lower = normal_cdf(lower, mean, sd);
    double normalization = cdf_upper - cdf_lower;

    // If normalization is practically zero, return 0 to match R behavior implicitly
    // (though R would technically divide by zero or epsilon)
    if (std::abs(normalization) < 1e-300) {
        return 0.0;
    }

    return normal_pdf(x, mean, sd) / normalization;
}

std::vector<SegmentCall> FastCall::assign_labels(
    const std::vector<std::vector<double>>& posteriors,
    const std::vector<double>& data
) {
    // Ported from LabelAss in LibraryFastCall.R
    
    static const std::array<int, N_STATES> CN_CALLS = {-2, -1, 0, 1, 2};
    static const std::array<int, N_STATES> ABSOLUTE_CN = {0, 1, 2, 3, 4};

    size_t n = posteriors.size();
    std::vector<SegmentCall> calls(n);

    for (size_t i = 0; i < n; ++i) {
        // R Code: indcall<-max.col(P0)
        // max.col with default ties.method="random" (usually first if stable impl)
        int max_state = 0;
        double max_prob = posteriors[i][0];

        for (int j = 1; j < N_STATES; ++j) {
            if (posteriors[i][j] > max_prob) {
                max_prob = posteriors[i][j];
                max_state = j;
            }
        }

        calls[i].cn_call = CN_CALLS[max_state];
        calls[i].absolute_cn = ABSOLUTE_CN[max_state];
        calls[i].probability = max_prob;
        calls[i].state_index = max_state;
        calls[i].segment_mean = data[i];
    }
    return calls;
}

double FastCall::compute_log_likelihood(
    const std::vector<double>& data,
    const std::vector<std::vector<double>>& posteriors
) {
    // Ported from EMFastCall in LibraryFastCall.R:
    // LikeliNew <- sum(PosteriorP(mdata,muvec,sdvec,prior) * prior)
    
    // IMPORTANT:
    // In R, 'PosteriorP(...)' returns an N x 5 matrix.
    // 'prior' is a vector of length 5.
    // The multiplication '*' in R performs element-wise multiplication with RECYCLING.
    // Since matrices in R are column-major, the 'prior' vector is recycled down the columns.
    // This means:
    //   Column 0 (State 0) is multiplied by prior[0], prior[1], prior[2]...
    //   Column 1 (State 1) continues the sequence.
    // This results in a mathematically non-standard weighting based on row index.
    // We strictly implement this R behavior here.

    double sum = 0.0;
    size_t n = data.size();

    // Iterate in Column-Major order to match R's recycling direction
    for (int j = 0; j < N_STATES; ++j) {
        for (size_t i = 0; i < n; ++i) {
            // Calculate flat index 'k' in R's column-major representation
            // k = j * n + i
            size_t flat_index = j * n + i;
            
            // Determine which prior index corresponds to this element
            size_t prior_index = flat_index % N_STATES;
            
            // Multiply posterior by the recycled prior
            sum += posteriors[i][j] * priors_[prior_index];
        }
    }
    
    return sum;
}

} // namespace fastcall
} // namespace excavator