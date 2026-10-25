// The loop shared by the kernels of the third and fourth derivatives of the
// expected information (src/d34e_<family>.cpp). A kernel returns, on the
// parameter scale, K components in the order of d3expected_names() or
// d4expected_names(); R sets the names (register_d34expected()).
//
// body(i, o) writes o[0], ..., o[K-1] for observation i; o writes straight
// into the output columns, so no buffer bounds K. The loop runs over
// d7::par_for(), and the result is identical at any thread count.

#ifndef DISTRIBUTIONS7_D34E_RUN_H
#define DISTRIBUTIONS7_D34E_RUN_H

#include <Rcpp.h>
#include <vector>
#include "d7_par.h"

namespace d34e {

// One row of the output: o[k] is component k at observation i.
struct Row {
    double* const* p;
    std::size_t i;
    double& operator[](int k) const { return p[k][i]; }
};

template <typename Body>
Rcpp::List run(int n, int threads, int threshold, int K, const Body& body) {
    std::vector<Rcpp::NumericVector> v(K);
    std::vector<double*> p(K);
    for (int k = 0; k < K; ++k) {
        v[k] = Rcpp::NumericVector(n);
        p[k] = v[k].begin();
    }
    double* const* pp = p.data();
    d7::par_for(n, threads, threshold, [&](std::size_t i) {
        body(i, Row{pp, i});
    });
    Rcpp::List out(K);
    for (int k = 0; k < K; ++k) out[k] = v[k];
    return out;
}

// A parameter read at observation i, recycled when it is a scalar.
struct Par {
    const double* x;
    bool s;
    explicit Par(const Rcpp::NumericVector& v) : x(v.begin()), s(v.size() == 1) {}
    double operator[](std::size_t i) const { return s ? x[0] : x[i]; }
};

}  // namespace d34e

#endif
