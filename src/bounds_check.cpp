#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Whether every entry of v satisfies is.finite(v) & v > lo & v < hi, the
// test of check_bounds_fast(), in one pass and without the four logical
// vectors the R expression allocates. NA and NaN fail is.finite() as they
// do in R. check_bounds_fast() calls the R expression only when this
// returns false, to build its message.
// [[Rcpp::export]]
bool bounds_ok_cpp(NumericVector v, double lo, double hi) {
    const double* x = v.begin();
    const R_xlen_t n = v.size();
    for (R_xlen_t i = 0; i < n; ++i) {
        const double a = x[i];
        if (!(std::isfinite(a) && a > lo && a < hi)) return false;
    }
    return true;
}
