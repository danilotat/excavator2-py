#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <htslib/sam.h>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

namespace py = pybind11;

namespace {
struct Boundary {
    std::int64_t position;
    std::size_t window;
    bool inclusive;
};

struct RegionPlan {
    std::string chromosome;
    std::size_t windows;
    std::int64_t begin, end;
    std::vector<Boundary> boundaries;

    RegionPlan(const std::string& chromosome_, const std::vector<std::int64_t>& starts,
               const std::vector<std::int64_t>& ends) : chromosome(chromosome_) {
        if (starts.empty() || starts.size() != ends.size())
            throw std::invalid_argument("counting requires matching nonempty window vectors");
        boundaries.reserve(2 * starts.size());
        for (std::size_t i = 0; i < starts.size(); ++i) {
            if (starts[i] > ends[i] || (i && starts[i] < starts[i - 1]))
                throw std::invalid_argument("counting requires windows ordered by start");
            boundaries.push_back({starts[i], i, false});
            boundaries.push_back({ends[i], i, true});
        }
        std::sort(boundaries.begin(), boundaries.end(), [](const Boundary& a, const Boundary& b) {
            return a.position < b.position ||
                   (a.position == b.position && a.inclusive < b.inclusive);
        });
        windows = starts.size();
        begin = std::max<std::int64_t>(1, starts.front()) - 1;
        end = std::max<std::int64_t>(1, *std::max_element(ends.begin(), ends.end()));
    }
};

class BamReader {
    std::string path;
    std::mutex mutex;
    std::unique_ptr<htsFile, decltype(&hts_close)> bam{nullptr, hts_close};
    std::unique_ptr<sam_hdr_t, decltype(&sam_hdr_destroy)> header{nullptr, sam_hdr_destroy};
    std::unique_ptr<hts_idx_t, decltype(&hts_idx_destroy)> index{nullptr, hts_idx_destroy};
    std::unique_ptr<bam1_t, decltype(&bam_destroy1)> read{nullptr, bam_destroy1};

public:
    explicit BamReader(const std::string& filename) : path(filename) {
        bam.reset(sam_open(path.c_str(), "rb"));
        if (!bam) throw std::invalid_argument("cannot open BAM: " + path);
        if (hts_get_format(bam.get())->format != ::bam)
            throw std::invalid_argument("preparation requires BAM; CRAM is not yet qualified");
        header.reset(sam_hdr_read(bam.get()));
        if (!header) throw std::invalid_argument("cannot read BAM header: " + path);
        index.reset(sam_index_load(bam.get(), path.c_str()));
        if (!index) throw std::invalid_argument("BAM preparation requires a readable index");
        read.reset(bam_init1());
        if (!read) throw std::bad_alloc();
    }

    void close() {
        std::lock_guard<std::mutex> lock(mutex);
        read.reset();
        index.reset();
        header.reset();
        bam.reset();
    }

    std::vector<std::int64_t> count(const RegionPlan& plan, int mapq) {
        std::lock_guard<std::mutex> lock(mutex);
        if (!bam) throw std::invalid_argument("BAM reader is closed");
        if (mapq < 0 || mapq > 255) throw std::invalid_argument("invalid MAPQ");
        const auto& chromosome = plan.chromosome;
        const auto& boundaries = plan.boundaries;
        int tid = sam_hdr_name2tid(header.get(), chromosome.c_str());
        if (tid < 0) throw std::invalid_argument("BAM lacks chromosome: " + chromosome);
        std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iterator(
            sam_itr_queryi(index.get(), tid, plan.begin, plan.end), hts_itr_destroy);
        if (!iterator) throw std::invalid_argument("cannot query chromosome: " + chromosome);
        std::vector<std::int64_t> result(plan.windows, 0);
        std::size_t next = 0;
        std::int64_t total = 0, previous = -1;
        auto record = [&](const Boundary& boundary) {
            result[boundary.window] += boundary.inclusive ? total : -total;
        };
        int status;
        while ((status = sam_itr_next(bam.get(), iterator.get(), read.get())) >= 0) {
            const auto& core = read->core;
            if ((core.flag & 1028) || core.qual < mapq) continue;
            const std::int64_t position = core.pos + 1;
            if (position < previous) throw std::invalid_argument("selected positions must be sorted");
            previous = position;
            // Prefix counts before starts and through ends preserve inclusive overlapping windows.
            while (next < boundaries.size() &&
                   (boundaries[next].position < position ||
                    (boundaries[next].position == position && !boundaries[next].inclusive))) {
                record(boundaries[next++]);
            }
            ++total;
        }
        if (status < -1) throw std::invalid_argument("error reading BAM: " + path);
        while (next < boundaries.size()) record(boundaries[next++]);
        return result;
    }
};
}  // namespace

void bind_reads(py::module_& module) {
    py::class_<RegionPlan>(module, "RegionPlan")
        .def(py::init<const std::string&, const std::vector<std::int64_t>&,
                      const std::vector<std::int64_t>&>(),
             py::arg("chromosome"), py::arg("starts"), py::arg("ends"),
             py::call_guard<py::gil_scoped_release>());
    py::class_<BamReader>(module, "BamReader")
        .def(py::init<const std::string&>(), py::arg("path"),
             py::call_guard<py::gil_scoped_release>())
        .def("count", &BamReader::count, py::arg("plan"), py::arg("mapq"),
             py::call_guard<py::gil_scoped_release>())
        .def("close", &BamReader::close, py::call_guard<py::gil_scoped_release>());
}
