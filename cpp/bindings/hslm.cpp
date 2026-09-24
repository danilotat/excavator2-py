#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <cmath>
#include <cstdint>
#include <limits>
#include "excavator2/hslm.hpp"

namespace py=pybind11;
using Array=py::array_t<double,py::array::c_style>;
namespace {
void vector(const Array& a, py::ssize_t size) {
    if (a.ndim()!=1 || a.size()!=size || reinterpret_cast<std::uintptr_t>(a.data())%alignof(double))
        throw py::value_error("HSLM requires aligned vectors with matching sizes");
    for (py::ssize_t i=0;i<size;++i) if (!std::isfinite(a.data()[i]))
        throw py::value_error("HSLM requires finite inputs");
}
py::tuple segment(const Array& values,const Array& means,double mi,double smu,double sepsilon,
                  const Array& eta,const Array& initial,bool trace) {
    const auto n=values.size(), k=means.size();
    if (n<1 || k<1 || k>std::numeric_limits<std::int32_t>::max())
        throw py::value_error("HSLM requires nonempty data and states");
    vector(values,n); vector(means,k); vector(eta,n-1); vector(initial,k);
    if (!std::isfinite(mi) || !std::isfinite(smu) || !std::isfinite(sepsilon) || smu<=0 || sepsilon<=0)
        throw py::value_error("HSLM deviations must be positive and finite");
    for (py::ssize_t i=0;i<n-1;++i) if (eta.data()[i]<=0 || eta.data()[i]>1)
        throw py::value_error("HSLM eta must be in (0,1]");
    const auto max=std::numeric_limits<py::ssize_t>::max()/sizeof(double);
    if (n>max/k || k>max/k || (trace && n-1>max/k/k))
        throw py::value_error("HSLM dimensions overflow buffer sizes");
    py::array_t<std::int32_t> path(n);
    Array transitions(trace ? std::vector<py::ssize_t>{n-1,k,k} : std::vector<py::ssize_t>{0});
    Array emissions(trace ? std::vector<py::ssize_t>{n,k} : std::vector<py::ssize_t>{0});
    Array scores(trace ? std::vector<py::ssize_t>{n,k} : std::vector<py::ssize_t>{0});
    py::array_t<std::int32_t> predecessors(trace ? std::vector<py::ssize_t>{n,k} : std::vector<py::ssize_t>{0});
    const auto* x=values.data(); const auto* mu=means.data(); const auto* e=eta.data();
    const auto* p=initial.data(); auto* out=path.mutable_data();
    auto* trans=trace ? transitions.mutable_data() : nullptr;
    auto* emit=trace ? emissions.mutable_data() : nullptr;
    auto* score=trace ? scores.mutable_data() : nullptr;
    auto* pred=trace ? predecessors.mutable_data() : nullptr;
    {
        py::gil_scoped_release release;
        excavator2::hslm::segment(n,k,x,mu,mi,smu,sepsilon,e,p,out,trans,emit,score,pred);
    }
    return py::make_tuple(path,transitions,emissions,scores,predecessors);
}
}
void bind_hslm(py::module_& module) {
    module.def("hslm_segment",&segment,py::arg("values").noconvert(),py::arg("means").noconvert(),
               py::arg("mi"),py::arg("smu"),py::arg("sepsilon"),py::arg("eta").noconvert(),
               py::arg("initial").noconvert(),py::arg("trace")=false);
}
