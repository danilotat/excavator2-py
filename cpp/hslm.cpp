#include "excavator2/hslm.hpp"
#include <cmath>
#include <vector>

namespace excavator2::hslm {
namespace {
// Literal ELNSUM from FastJointSLMLibraryI.f; retain its evaluation order.
double elnsum(double x, double y) {
    return x > y ? x + std::log(1 + std::exp(y-x)) : y + std::log(1 + std::exp(x-y));
}
}
void segment(std::size_t n, std::size_t k, const double* values, const double* means,
             double mi, double smu, double sepsilon, const double* eta,
             const double* initial, std::int32_t* path, double* transitions,
             double* emissions, double* scores, std::int32_t* predecessors) {
    std::vector<double> g(k), previous(k), current(k), block(k*k);
    std::vector<std::int32_t> history(n*k, 0);
    for (std::size_t j=0; j<k; ++j) {
        const double difference = means[j]-mi;
        g[j] = -(difference*difference/(2*(smu*smu)));
    }
    double norm=g[0];
    for (std::size_t j=1; j<k; ++j) norm=elnsum(norm,g[j]);
    for (auto& value:g) value-=norm;
    // The unsuffixed Fortran PI literal is rounded to REAL*4, then promoted.
    const double constant=std::log(1/(std::sqrt(2*3.1415927410125732421875)*sepsilon));
    for (std::size_t t=0; t<n; ++t) {
        if (t>0) {
            const double jump=std::log(eta[t-1]), stay=std::log(1-eta[t-1]);
            for (std::size_t source=0; source<k; ++source) {
                for (std::size_t dest=0; dest<k; ++dest) {
                    const double change=jump+g[dest];
                    block[source*k+dest]=source==dest ? elnsum(stay,change) : change;
                    if (transitions) transitions[(t-1)*k*k+source*k+dest]=block[source*k+dest];
                }
            }
        }
        for (std::size_t j=0; j<k; ++j) {
            const double z=(values[t]-means[j])/sepsilon;
            const double emission=constant+(-0.5*(z*z));
            if (emissions) emissions[t*k+j]=emission;
            if (t==0) current[j]=initial[j]+emission;
            else {
                double best=previous[0]+block[j];
                std::int32_t index=1;
                for (std::size_t source=1; source<k; ++source) {
                    const double candidate=previous[source]+block[source*k+j];
                    if (candidate>best) { best=candidate; index=static_cast<std::int32_t>(source+1); }
                }
                history[t*k+j]=index;
                current[j]=best+emission;
            }
            if (scores) scores[t*k+j]=current[j];
            if (predecessors) predecessors[t*k+j]=history[t*k+j];
        }
        previous.swap(current);
    }
    std::size_t best=0;
    for (std::size_t j=1; j<k; ++j) if (previous[j]>previous[best]) best=j;
    path[n-1]=static_cast<std::int32_t>(best+1);
    for (std::size_t t=n-1; t>0; --t) path[t-1]=history[t*k+path[t]-1];
}
}
