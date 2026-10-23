#include <Rcpp.h>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// Fifth-order derivatives of the BB2 log-mass in the shapes (alpha, beta),
// size n, written out. With s = alpha + beta the log-mass is
//   l = lgamma(y+alpha) + lgamma(n-y+beta) - lgamma(n+s)
//       - lgamma(alpha) - lgamma(beta) + lgamma(s) + const.
// Every derivative of order five is a polygamma difference of order four
// (derive.py, sympy, checks the pairing):
//   pure alpha:  psi4(y+alpha) - psi4(alpha) + [psi4(s) - psi4(n+s)]
//   mixed:       psi4(s) - psi4(n+s)
//   pure beta:   psi4(n-y+beta) - psi4(beta) + [psi4(s) - psi4(n+s)]
// The component with b derivatives in beta and 5 - b in alpha is named
// deriv_names(c("alpha", "beta"), 5) in that order.
//
// The cut is the one deriv4 uses, s = alpha + beta > 10 (bb_direct_cut() of
// R/betabinom2_distrib.R). Above it the polygamma differences are replaced
// by the sums over the support that bb2_direct() of src/pt_betabinom.h writes
// for orders up to four, here at order five:
//   c_5 = +24,
//   pure alpha:  24 [ sum_{i<y} 1/(alpha+i)^5 - sum_{i<n} 1/(s+i)^5 ],
//   pure beta:   24 [ sum_{i<n-y} 1/(beta+i)^5 - sum_{i<n} 1/(s+i)^5 ],
//   mixed:       -24 sum_{i<n} 1/(s+i)^5.
// Each sum has terms of one sign, and the sums cost n terms an observation,
// as the deriv4 sums do. Below the cut the polygamma form costs O(1), and
// every argument there is at most n + 10.

static inline void d5bb2_row(double y, double n, double al, double be, double* out) {
    const double s = al + be;
    if (s > 10.0) {
        double SA = 0.0, SB = 0.0, ST = 0.0;
        const long ky = (long) y, kv = (long) (n - y), kn = (long) n;
        for (long i = 0; i < kn; ++i) {
            const double di = (double) i;
            const double it = 1.0 / (s + di);
            const double it2 = it * it;
            ST += it2 * it2 * it;                         // 1/(s+i)^5
            if (i < ky) {
                const double ia = 1.0 / (al + di);
                const double ia2 = ia * ia;
                SA += ia2 * ia2 * ia;                     // 1/(alpha+i)^5
            }
            if (i < kv) {
                const double ib = 1.0 / (be + di);
                const double ib2 = ib * ib;
                SB += ib2 * ib2 * ib;                     // 1/(beta+i)^5
            }
        }
        for (int b = 0; b <= 5; ++b) {
            double v = -ST;
            if (b == 0) v += SA;
            if (b == 5) v += SB;
            out[b] = 24.0 * v;
        }
        return;
    }
    // psi^(4)(s) - psi^(4)(n+s), the part every component carries
    const double mixed = R::psigamma(s, 4.0) - R::psigamma(n + s, 4.0);
    out[0] = R::psigamma(y + al, 4.0) - R::psigamma(al, 4.0) + mixed;
    out[1] = mixed;
    out[2] = mixed;
    out[3] = mixed;
    out[4] = mixed;
    out[5] = R::psigamma(n - y + be, 4.0) - R::psigamma(be, 4.0) + mixed;
}

// [[Rcpp::export]]
List betabinom2_deriv5_cpp(NumericVector y, NumericVector alpha, NumericVector beta,
                           double size) {
    const int ny = y.size(), na = alpha.size(), nb = beta.size();
    const int n = std::max(ny, std::max(na, nb));
    NumericVector c_aaaaa(n), c_aaaab(n), c_aaabb(n), c_aabbb(n), c_abbbb(n), c_bbbbb(n);

    for (int i = 0; i < n; ++i) {
        double out[6];
        d5bb2_row(y[i % ny], size, alpha[i % na], beta[i % nb], out);
        c_aaaaa[i] = out[0];
        c_aaaab[i] = out[1];
        c_aaabb[i] = out[2];
        c_aabbb[i] = out[3];
        c_abbbb[i] = out[4];
        c_bbbbb[i] = out[5];
    }

    return List::create(
        Named("alpha_alpha_alpha_alpha_alpha") = c_aaaaa,
        Named("alpha_alpha_alpha_alpha_beta") = c_aaaab,
        Named("alpha_alpha_alpha_beta_beta") = c_aaabb,
        Named("alpha_alpha_beta_beta_beta") = c_aabbb,
        Named("alpha_beta_beta_beta_beta") = c_abbbb,
        Named("beta_beta_beta_beta_beta") = c_bbbbb
    );
}
