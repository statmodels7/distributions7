#include <Rcpp.h>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// Fifth-order derivatives of the BB1 log-mass in (mu, sigma), size n, written
// out. The shapes are alpha = mu/sigma and beta = (1 - mu)/sigma, with
// s = alpha + beta = 1/sigma, and the log-mass is the BB2 one at those shapes
// (betabinom2_deriv5_cpp() of d5_betabinom2.cpp). The derivative in (mu,
// sigma) follows from the chain of the map, which derive.py (sympy) carries out:
// each polygamma pair enters with the coefficient it has in the map's chain,
// and the pairs are
//   Dal_k = psi^(k)(y+alpha) - psi^(k)(alpha),
//   Dbe_k = psi^(k)(n-y+beta) - psi^(k)(beta),
//   Ds_k  = psi^(k)(n+s) - psi^(k)(s), s = 1/sigma.
// The coefficients are powers of sigma^-1 and polynomials in mu, and written
// below in the factored form derive.py printed. The component with b
// derivatives in sigma and 5 - b in mu is the b-th in the order of
// deriv_names(c("mu", "sigma"), 5).
//
// The pair differences psi^(k)(x+m) - psi^(k)(x) are sums over the shift (d5bb1_pair)
// when m <= 32, which is every shift a small count gives, and polygamma differences
// above that.
//
// The cut is sigma < 0.1, bb1_direct_cut() of R/betabinom2_distrib.R. Below it
// the map's own derivatives cancel among themselves as sigma goes to zero, so
// the components are the sums over the support of pt_betabinom.h's bb1_direct,
// written here at order five:
//   c_5 = +24, and with SU_b = sum_{i<y} i^b/(mu+i sigma)^5,
//   SV_b = sum_{i<n-y} i^b/(1-mu+i sigma)^5, SW = sum_{i<n} i^5/(1+i sigma)^5,
//   component b (a = 5 - b mu-derivatives):
//     24 [ SU_b + (-1)^a SV_b - [a = 0] SW ].
// Above it the polygamma pairs cost O(1), and every argument is at most
// n + 1/sigma <= n + 10, since mu lies in (0, 1).

static inline void d5bb1_direct(double y, double n, double mu, double s, double* out) {
    double SU[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double SV[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double SW = 0.0;
    const long ky = (long) y, kv = (long) (n - y), kn = (long) n;
    const double nu = 1.0 - mu;
    for (long i = 0; i < kn; ++i) {
        const double di = (double) i;
        if (i < ky) {
            const double iu = 1.0 / (mu + di * s);
            const double iu2 = iu * iu;
            const double base = iu2 * iu2 * iu;                  // 1/u^5
            double pw = 1.0;
            for (int b = 0; b <= 5; ++b) { SU[b] += pw * base; pw *= di; }
        }
        if (i < kv) {
            const double iv = 1.0 / (nu + di * s);
            const double iv2 = iv * iv;
            const double base = iv2 * iv2 * iv;                  // 1/v^5
            double pw = 1.0;
            for (int b = 0; b <= 5; ++b) { SV[b] += pw * base; pw *= di; }
        }
        const double iw = 1.0 / (1.0 + di * s);
        const double iw2 = iw * iw;
        SW += di * di * di * di * di * iw2 * iw2 * iw;           // i^5 / (1+i s)^5
    }
    for (int b = 0; b <= 5; ++b) {
        const int a = 5 - b;
        double v = SU[b] + ((a % 2 == 0) ? SV[b] : -SV[b]);
        if (a == 0) v -= SW;
        out[b] = 24.0 * v;
    }
}

// psi^(k)(x + m) - psi^(k)(x) for an integer shift m. For m <= 32 it is the
// sum of (-1)^k k! / (x+i)^(k+1) over i < m, whose terms all have one sign, so
// nothing cancels; above that it is the polygamma difference. The sum is the
// device psi_diff1() of pt_betabinom.h uses at order one.
static inline double d5bb1_pair(int k, double m, double x) {
    if (m >= 0.0 && m <= 32.0 && m == std::floor(m)) {
        double acc = 0.0;
        for (int i = 0; i < (int) m; ++i) {
            const double t = 1.0 / (x + (double) i);
            double p = t;
            for (int j = 0; j < k; ++j) p *= t;             // t^(k+1)
            acc += p;
        }
        double fk = 1.0;
        for (int j = 2; j <= k; ++j) fk *= (double) j;      // k!
        return ((k % 2 == 0) ? 1.0 : -1.0) * fk * acc;
    }
    return R::psigamma(x + m, (double) k) - R::psigamma(x, (double) k);
}

static inline void d5bb1_row(double y, double n, double mu, double sg, double* out) {
    if (sg < 0.1) {
        d5bb1_direct(y, n, mu, sg, out);
        return;
    }
    const double al = mu / sg, be = (1.0 - mu) / sg, s = 1.0 / sg;
    // the pair differences, k = 0..4
    double Da[5], Db[5], Ds[5];
    for (int k = 0; k < 5; ++k) {
        Da[k] = d5bb1_pair(k, y, al);
        Db[k] = d5bb1_pair(k, n - y, be);
        Ds[k] = d5bb1_pair(k, n, s);
    }
    const double s2 = sg * sg;
    const double s5 = s2 * s2 * sg;
    const double s6 = s5 * sg, s7 = s6 * sg, s8 = s7 * sg, s9 = s8 * sg, s10 = s9 * sg;
    const double m1 = mu - 1.0;
    const double mu2 = mu * mu, mu3 = mu2 * mu, mu4 = mu3 * mu, mu5 = mu4 * mu;
    const double m12 = m1 * m1, m13 = m12 * m1, m14 = m13 * m1, m15 = m14 * m1;

    // mu^5
    out[0] = (Da[4] - Db[4]) / s5;
    // mu^4 sigma
    out[1] = (-4.0 * Da[3] - 4.0 * Db[3]) / s5
             + (-mu * Da[4] + m1 * Db[4]) / s6;
    // mu^3 sigma^2
    out[2] = (12.0 * Da[2] - 12.0 * Db[2]) / s5
             + (8.0 * mu * Da[3] + 8.0 * m1 * Db[3]) / s6
             + (mu2 * Da[4] - m12 * Db[4]) / s7;
    // mu^2 sigma^3
    out[3] = (-24.0 * Da[1] - 24.0 * Db[1]) / s5
             + (-36.0 * mu * Da[2] + 36.0 * m1 * Db[2]) / s6
             + (-12.0 * mu2 * Da[3] - 12.0 * m12 * Db[3]) / s7
             + (-mu3 * Da[4] + m13 * Db[4]) / s8;
    // mu sigma^4
    out[4] = (24.0 * Da[0] - 24.0 * Db[0]) / s5
             + (96.0 * mu * Da[1] + 96.0 * m1 * Db[1]) / s6
             + (72.0 * mu2 * Da[2] - 72.0 * m12 * Db[2]) / s7
             + (16.0 * mu3 * Da[3] + 16.0 * m13 * Db[3]) / s8
             + (mu4 * Da[4] - m14 * Db[4]) / s9;
    // sigma^5
    out[5] = (-120.0 * mu * Da[0] + 120.0 * m1 * Db[0]) / s6
             + (-240.0 * mu2 * Da[1] - 240.0 * m12 * Db[1]) / s7
             + (-120.0 * mu3 * Da[2] + 120.0 * m13 * Db[2]) / s8
             + (-20.0 * mu4 * Da[3] - 20.0 * m14 * Db[3]) / s9
             + (-mu5 * Da[4] + m15 * Db[4]) / s10
             + (120.0 * Ds[0] / s6 + 240.0 * Ds[1] / s7 + 120.0 * Ds[2] / s8
                + 20.0 * Ds[3] / s9 + Ds[4] / s10);
}

// [[Rcpp::export]]
List betabinom1_deriv5_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                           double size) {
    const int ny = y.size(), nm = mu.size(), ns = sigma.size();
    const int n = std::max(ny, std::max(nm, ns));
    NumericVector c_mmmmm(n), c_mmmms(n), c_mmmss(n), c_mmsss(n), c_msssss(n), c_sssss(n);

    for (int i = 0; i < n; ++i) {
        double out[6];
        d5bb1_row(y[i % ny], size, mu[i % nm], sigma[i % ns], out);
        c_mmmmm[i] = out[0];
        c_mmmms[i] = out[1];
        c_mmmss[i] = out[2];
        c_mmsss[i] = out[3];
        c_msssss[i] = out[4];
        c_sssss[i] = out[5];
    }

    return List::create(
        Named("mu_mu_mu_mu_mu") = c_mmmmm,
        Named("mu_mu_mu_mu_sigma") = c_mmmms,
        Named("mu_mu_mu_sigma_sigma") = c_mmmss,
        Named("mu_mu_sigma_sigma_sigma") = c_mmsss,
        Named("mu_sigma_sigma_sigma_sigma") = c_msssss,
        Named("sigma_sigma_sigma_sigma_sigma") = c_sssss
    );
}
