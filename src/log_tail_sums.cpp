#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// The partition sums of quad_log_tail_derivs(), one row of quantiles at a
// time. For every component c (a multi-index of order 1 to 4 in canonical
// form), the complete Bell polynomial in the partials of D at a node obeys
//   B(c) = sum over the subsets s of c's positions that contain the first
//          of D(s) B(c \ s),        B(empty) = 1,
// and the cumulants of the weighted means M(c) = sum_j w_j B_j(c) obey
//   M(c) = sum over the same subsets of kappa(s) M(c \ s).
// The pairs (s, c \ s) are enumerated once in R: for component c they are
// ps[k], pr[k] for k in [pstart[c], pstart[c + 1]), each the index of a
// component, pr = -1 for the empty rest. Components are ordered so that
// every subset comes before the component it is part of.

// tabs holds one vector per component, the log-density partial at the n
// quantiles followed by the n x K nodes (column-major), and D at row i and
// node j is tabs[c][n + i + n j] - tabs[c][i].

// [[Rcpp::export]]
NumericMatrix log_tail_cumulants_cpp(List tabs, NumericMatrix W,
                                     IntegerVector pstart, IntegerVector ps,
                                     IntegerVector pr) {
    const int n = W.nrow(), K = W.ncol();
    const int C = pstart.size() - 1;
    std::vector<const double*> tv(C);
    for (int c = 0; c < C; ++c) tv[c] = REAL(tabs[c]);
    NumericMatrix out(n, C);
    // the means M(c) for every row, accumulated node by node so that the
    // quantiles of one node, contiguous in tabs and W, are read in order
    std::vector<double> M((std::size_t) n * C, 0.0), B(C), D(C), kap(C);
    for (int j = 0; j < K; ++j) {
        const double* wj = &W(0, j);
        for (int i = 0; i < n; ++i) {
            const double w = wj[i];
            if (w == 0.0) continue;
            const std::size_t at = (std::size_t) n + i + (std::size_t) n * j;
            for (int c = 0; c < C; ++c) D[c] = tv[c][at] - tv[c][i];
            double* Mi = &M[(std::size_t) i * C];
            for (int c = 0; c < C; ++c) {
                double s = 0.0;
                for (int k = pstart[c]; k < pstart[c + 1]; ++k) {
                    const double d = D[ps[k]];
                    s += (pr[k] < 0) ? d : d * B[pr[k]];
                }
                B[c] = s;
                Mi[c] += w * s;
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        const double* Mi = &M[(std::size_t) i * C];
        for (int c = 0; c < C; ++c) {
            double s = Mi[c];
            for (int k = pstart[c]; k < pstart[c + 1]; ++k) {
                if (pr[k] < 0) continue;
                s -= kap[ps[k]] * Mi[pr[k]];
            }
            kap[c] = s;
            out(i, c) = s;
        }
    }
    return out;
}
