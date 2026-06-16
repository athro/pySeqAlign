// Fast C++ Needleman-Wunsch affine-gap aligner for pyREAL boosting.
//
// Mirrors pyreal.core.nw_affine.NeedlemanWunschAffine EXACTLY (3-matrix M/Ix/Iy
// DP, same border init, same argmax tie-breaking, same traceback + gap counting)
// so the Python and C++ backends are interchangeable. Adapted from the pyAlign2
// AlignerAffine structure, but YAP-free: scores come from a flat distance matrix
// (the boosting reward matrix), set once per reward update; align() is then a
// pure numeric kernel callable thousands of times without recomputing scores.
//
// Build: see build_cpp_aligner.sh
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <limits>
#include <algorithm>
#include <stdexcept>

namespace py = pybind11;

static const int M = 0, IX = 1, IY = 2;

struct AlignResult {
    double score = 0.0;
    std::vector<int> query;   // aligned query (0 = gap)
    std::vector<int> target;  // aligned target (0 = gap)
    int gap_opens = 0;
    int gap_extensions = 0;
    int length = 0;
};

class CppAligner {
public:
    // num_ids = number of distinct atom ids (gap id 0 excluded). The score
    // matrix is (num_ids+1) x (num_ids+1), row-major, indexed by atom id.
    CppAligner(int num_ids, double gap_open, double gap_extend)
        : n1_(num_ids + 1), gap_open_(gap_open), gap_extend_(gap_extend),
          mat_((size_t)(num_ids + 1) * (num_ids + 1), 0.0) {}

    void set_matrix(const std::vector<double>& flat) {
        if (flat.size() != mat_.size())
            throw std::runtime_error("score matrix size mismatch");
        mat_ = flat;
    }

    inline double score(int a, int b) const { return mat_[(size_t)a * n1_ + b]; }

    AlignResult align(const std::vector<int>& q, const std::vector<int>& t) const {
        const int n = (int)q.size(), m = (int)t.size();
        const double NEG = -std::numeric_limits<double>::infinity();
        const double d = gap_open_, e = gap_extend_;

        // F[k][i][j], B stores (from_k, from_i, from_j)
        std::vector<double> F((size_t)3 * (n + 1) * (m + 1), NEG);
        std::vector<int> B((size_t)3 * (n + 1) * (m + 1) * 3, -1);
        auto Fi = [&](int k, int i, int j) -> double& {
            return F[((size_t)k * (n + 1) + i) * (m + 1) + j];
        };
        auto Bi = [&](int k, int i, int j, int c) -> int& {
            return B[(((size_t)k * (n + 1) + i) * (m + 1) + j) * 3 + c];
        };
        auto setB = [&](int k, int i, int j, int fk, int fi, int fj) {
            Bi(k, i, j, 0) = fk; Bi(k, i, j, 1) = fi; Bi(k, i, j, 2) = fj;
        };

        Fi(M, 0, 0) = 0.0;
        for (int i = 1; i <= n; ++i) {
            Fi(IX, i, 0) = (i > 1) ? Fi(IX, i - 1, 0) + e : d;
            setB(IX, i, 0, IX, i - 1, 0);
        }
        for (int j = 1; j <= m; ++j) {
            Fi(IY, 0, j) = (j > 1) ? Fi(IY, 0, j - 1) + e : d;
            setB(IY, 0, j, IY, 0, j - 1);
        }

        for (int i = 1; i <= n; ++i) {
            for (int j = 1; j <= m; ++j) {
                double s = score(q[i - 1], t[j - 1]);
                double cm[3] = {Fi(M, i - 1, j - 1) + s, Fi(IX, i - 1, j - 1) + s, Fi(IY, i - 1, j - 1) + s};
                int bk = argmax3(cm);
                Fi(M, i, j) = cm[bk]; setB(M, i, j, bk, i - 1, j - 1);

                double cx[3] = {Fi(M, i - 1, j) + d, Fi(IX, i - 1, j) + e, Fi(IY, i - 1, j) + d};
                bk = argmax3(cx);
                Fi(IX, i, j) = cx[bk]; setB(IX, i, j, bk, i - 1, j);

                double cy[3] = {Fi(M, i, j - 1) + d, Fi(IY, i, j - 1) + e, Fi(IX, i, j - 1) + d};
                bk = argmax3(cy);
                Fi(IY, i, j) = cy[bk]; setB(IY, i, j, bk, i, j - 1);
            }
        }

        double ends[3] = {Fi(M, n, m), Fi(IX, n, m), Fi(IY, n, m)};
        int best = argmax3(ends);

        AlignResult r;
        r.score = ends[best];
        int k = best, i = n, j = m, prev_k = -1;
        while (i > 0 || j > 0) {
            int fk = Bi(k, i, j, 0), fi = Bi(k, i, j, 1), fj = Bi(k, i, j, 2);
            if (fi < 0 || fj < 0) break;
            if (k == M) { r.query.push_back(q[i - 1]); r.target.push_back(t[j - 1]); }
            else if (k == IX) {
                r.query.push_back(q[i - 1]); r.target.push_back(0);
                if (prev_k != IX) r.gap_opens++; else r.gap_extensions++;
            } else {
                r.query.push_back(0); r.target.push_back(t[j - 1]);
                if (prev_k != IY) r.gap_opens++; else r.gap_extensions++;
            }
            prev_k = k; k = fk; i = fi; j = fj;
        }
        std::reverse(r.query.begin(), r.query.end());
        std::reverse(r.target.begin(), r.target.end());
        r.length = (int)r.query.size();
        return r;
    }

private:
    int n1_;
    double gap_open_, gap_extend_;
    std::vector<double> mat_;

    static inline int argmax3(const double v[3]) {
        if (v[0] >= v[1]) return (v[0] >= v[2]) ? 0 : 2;
        return (v[1] >= v[2]) ? 1 : 2;
    }
};

PYBIND11_MODULE(cpp_aligner, mod) {
    py::class_<AlignResult>(mod, "AlignResult")
        .def_readonly("score", &AlignResult::score)
        .def_readonly("query", &AlignResult::query)
        .def_readonly("target", &AlignResult::target)
        .def_readonly("gap_opens", &AlignResult::gap_opens)
        .def_readonly("gap_extensions", &AlignResult::gap_extensions)
        .def_readonly("length", &AlignResult::length);
    py::class_<CppAligner>(mod, "CppAligner")
        .def(py::init<int, double, double>())
        .def("set_matrix", &CppAligner::set_matrix)
        .def("align", &CppAligner::align);
}
