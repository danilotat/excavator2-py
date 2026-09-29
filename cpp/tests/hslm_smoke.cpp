#include "excavator2/hslm.hpp"
#include <cassert>
#include <cmath>
#include <vector>

int main() {
    for (std::size_t n : {1, 17, 10003}) {
        const std::size_t k=21;
        std::vector<double> values(n), means(k), eta(n-1,0.0001), initial(k,std::log(1.0/k));
        for (std::size_t j=0;j<k;++j) means[j]=-1+0.1*j;
        for (std::size_t i=0;i<n;++i) values[i]=std::sin(i*0.1);
        std::vector<std::int32_t> path(n), traced(n), predecessors(n*k);
        std::vector<double> transitions((n-1)*k*k), emissions(n*k), scores(n*k);
        excavator2::hslm::segment(n,k,values.data(),means.data(),0,0.1,0.3,eta.data(),
                                 initial.data(),path.data(),nullptr,nullptr,nullptr,nullptr);
        excavator2::hslm::segment(n,k,values.data(),means.data(),0,0.1,0.3,eta.data(),
                                 initial.data(),traced.data(),transitions.data(),emissions.data(),
                                 scores.data(),predecessors.data());
        assert(path==traced);
        for (auto state:path) assert(state>=1 && state<=21);
        for (auto score:scores) assert(std::isfinite(score));
    }
}
