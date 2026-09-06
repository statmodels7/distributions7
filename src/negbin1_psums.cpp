#include <Rcpp.h>
#include "d7_par.h"
using namespace Rcpp;

// THE POWER SUMS THE NB1 HIGHER DERIVATIVES ARE WRITTEN IN.
//
// In the size r = mu/theta the NB1 log-mass is G(r) + r B(theta) + C(theta),
// and every derivative of G is a polygamma differenced at the shift y.  As
// theta goes to zero the family tends to the Poisson, r runs away, and the
// assembly divides by theta^(a+b): the terms r^j G^(a+j)(r) are then all of
// one size and sum to something far smaller, so the answer is lost among them.
// Measured at mu = 4, y = 3, the third derivative in theta reads 2.770e-01 at
// theta = 1e-3 where it is 2.797e-01, -2.37 at 1e-4, -1.97e+09 at 1e-6 and
// -1.46e+17 at 1e-8; the fourth reads 8.97 at 1e-3 where it is -1.59.
//
// The size need not appear at all.  Since lgamma(y+r) - lgamma(r) is
// sum_{i<y} log(r + i) and r = mu/theta, the y log(theta) that each such term
// carries cancels EXACTLY against the C(theta) part, leaving
//
//   l = sum_{i<y} log(mu + i theta) - log(y!)
//       - (mu/theta) log1p(theta) - y log1p(theta),
//
// whose natural variable is mu + i theta and does not run away.  Verified
// against dnbinom to between 0 and 2.2e-11 over counts, means and dispersions
// from 1e-6 to 20, the last being R's own loss of the Poisson limit.
//
// Every derivative of the first term is then
//
//   d^a/dmu^a d^b/dtheta^b sum_{i<y} log(mu + i theta)
//     = (-1)^(a+b-1) (a+b-1)! sum_{i<y} i^b / (mu + i theta)^(a+b),
//
// so ONE pass over i serves every component of an order: the denominator power
// is a + b, which is the order, and only the numerator power b changes.  This
// returns those sums, column b holding sum_{i<y} i^b/(mu + i theta)^order for
// b = 0..order.
//
// The cost is y terms an observation against the O(1) of the polygamma route,
// which is why the caller takes this branch only below a measured crossover in
// theta.  On a real design of 57600 counts summing to 6.24e6 the whole pass is
// 6.24e6 terms, a few milliseconds.
//
// i^b is exact in a double for every count a design carries: at b = 4 it
// passes 2^53 only above i = 9.7e3, and above that the loop is already the
// dominant cost rather than the accuracy.

// [[Rcpp::export]]
NumericMatrix negbin1_psums_cpp(NumericVector y, NumericVector mu,
                                NumericVector theta, int order,
                                int threads = 1) {
    const int n = y.size();
    const int nb = order + 1;
    NumericMatrix out(n, nb);
    const bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t idx) {
        const int i = (int) idx;
        const double m = mu_s ? mu[0] : mu[i];
        const double th = th_s ? theta[0] : theta[i];
        const double yi = y[i];
        double acc[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        double t = m;                       // mu + k theta at k = 0
        for (double k = 0.0; k < yi; k += 1.0) {
            double d = t;
            for (int p = 1; p < order; ++p) d *= t;
            const double inv = 1.0 / d;
            double kb = 1.0;                // k^b
            for (int b = 0; b < nb; ++b) {
                acc[b] += kb * inv;
                kb *= k;
            }
            t += th;
        }
        for (int b = 0; b < nb; ++b) out(i, b) = acc[b];
    });
    return out;
}
