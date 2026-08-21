#include <Rcpp.h>
#include "d7_par.h"

using namespace Rcpp;

// The near-zero branch of the generalized Pareto's third and fourth
// derivatives. gpd_components() writes each component there as a series in
// u = xi z, and with the two elementwise powers the loop used to raise
// collapsed into one -- xi^(k-b) z^(k+1) = u^(k-b) z^(b+1), so z^(b+1)/sigma^a
// leaves the loop -- what remains is a POLYNOMIAL IN u with scalar
// coefficients, one per element. Forty-one terms over the elements below the
// cut, five components at fourth order, is a scalar recursion of a few
// million steps: measured, the R form costs 510 ms at n = 20000 and is 97 per
// cent of the derivative, where the branch it complements costs 12.
//
// The evaluation is Horner from the highest power down, which for a series
// whose terms decay adds the smallest first and is the better summation order
// as well as the shorter one; the R loop summed largest-first.
//
// The body is arithmetic and nothing else -- no Rmath, no allocation -- so it
// is admissible in a worker, and the decomposition is over the elements of
// the output.
// [[Rcpp::export]]
NumericVector gpd_poly_cpp(NumericVector u, NumericVector coef, int threads) {
    const R_xlen_t n = u.size();
    const int m = (int) coef.size();
    NumericVector out(n);
    const double* up = u.begin();
    const double* cp = coef.begin();
    double* op = out.begin();
    if (m == 0) return out;
    d7::par_for((int) n, threads, d7::kMinCostly, [&](std::size_t i) {
        const double x = up[i];
        double acc = cp[m - 1];
        for (int j = m - 2; j >= 0; --j) acc = acc * x + cp[j];
        op[i] = acc;
    });
    return out;
}
