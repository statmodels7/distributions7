#include <Rcpp.h>
using namespace Rcpp;

// Beta-binomial with size n, written in the mean proportion mu and the
// dispersion sigma. The shape parameters are A = mu/sigma and
// B = (1-mu)/sigma, so S = A + B = 1/sigma, and
//   l = lchoose(n, y) + lgamma(y+A) + lgamma(n-y+B) - lgamma(n+S)
//                     - lgamma(A) - lgamma(B) + lgamma(S).
// The whole dependence on the parameters is through A and B, where every
// derivative is a difference of polygammas:
//   l_A = psi(y+A) - psi(A) - [psi(n+S) - psi(S)],
// and likewise for B, with l_AB carrying only the S part. The chain rule to
// (mu, sigma) is then two variables deep and written out below.
//
// Var(Y) = n mu (1-mu) [1 + (n-1) sigma/(1+sigma)], so sigma -> 0 recovers
// the binomial and the family is overdispersed for every sigma > 0.

struct BBderiv {
    double lA, lB;              // first order in the shapes
    double lAA, lAB, lBB;       // second order
};

static inline BBderiv bb_shape_derivs(double y, double n, double A, double B) {
    double S = A + B;
    BBderiv d;
    double dS  = R::digamma(n + S) - R::digamma(S);
    double dS1 = R::trigamma(n + S) - R::trigamma(S);
    d.lA  = R::digamma(y + A) - R::digamma(A) - dS;
    d.lB  = R::digamma(n - y + B) - R::digamma(B) - dS;
    d.lAA = R::trigamma(y + A) - R::trigamma(A) - dS1;
    d.lBB = R::trigamma(n - y + B) - R::trigamma(B) - dS1;
    d.lAB = -dS1;
    return d;
}

// The shapes and their first two derivatives in (mu, sigma).
struct BBmap {
    double A, B;
    double Am, As, Amm, Ams, Ass;
    double Bm, Bs, Bmm, Bms, Bss;
};

static inline BBmap bb_map(double mu, double s) {
    BBmap m;
    double s2 = s * s, s3 = s2 * s;
    m.A = mu / s;          m.B = (1.0 - mu) / s;
    m.Am = 1.0 / s;        m.Bm = -1.0 / s;
    m.As = -mu / s2;       m.Bs = -(1.0 - mu) / s2;
    m.Amm = 0.0;           m.Bmm = 0.0;
    m.Ams = -1.0 / s2;     m.Bms = 1.0 / s2;
    m.Ass = 2.0 * mu / s3; m.Bss = 2.0 * (1.0 - mu) / s3;
    return m;
}

// [[Rcpp::export]]
List betabinom_gradient_cpp(NumericVector y, NumericVector mu,
                            NumericVector sigma, double size) {
    int n = y.size();
    NumericVector g_mu(n), g_sigma(n);
    bool mu_s = (mu.size() == 1), si_s = (sigma.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = mu_s ? mu[0] : mu[i];
        double s = si_s ? sigma[0] : sigma[i];
        BBmap mp = bb_map(m, s);
        BBderiv d = bb_shape_derivs(y[i], size, mp.A, mp.B);
        g_mu[i]    = d.lA * mp.Am + d.lB * mp.Bm;
        g_sigma[i] = d.lA * mp.As + d.lB * mp.Bs;
    }
    return List::create(Named("mu") = g_mu, Named("sigma") = g_sigma);
}

// [[Rcpp::export]]
List betabinom_hessian_cpp(NumericVector y, NumericVector mu,
                           NumericVector sigma, double size) {
    int n = y.size();
    NumericVector h_mm(n), h_ms(n), h_ss(n);
    bool mu_s = (mu.size() == 1), si_s = (sigma.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = mu_s ? mu[0] : mu[i];
        double s = si_s ? sigma[0] : sigma[i];
        BBmap mp = bb_map(m, s);
        BBderiv d = bb_shape_derivs(y[i], size, mp.A, mp.B);

        h_mm[i] = d.lAA * mp.Am * mp.Am
                + 2.0 * d.lAB * mp.Am * mp.Bm
                + d.lBB * mp.Bm * mp.Bm
                + d.lA * mp.Amm + d.lB * mp.Bmm;
        h_ms[i] = d.lAA * mp.Am * mp.As
                + d.lAB * (mp.Am * mp.Bs + mp.As * mp.Bm)
                + d.lBB * mp.Bm * mp.Bs
                + d.lA * mp.Ams + d.lB * mp.Bms;
        h_ss[i] = d.lAA * mp.As * mp.As
                + 2.0 * d.lAB * mp.As * mp.Bs
                + d.lBB * mp.Bs * mp.Bs
                + d.lA * mp.Ass + d.lB * mp.Bss;
    }
    return List::create(Named("mu_mu") = h_mm,
                        Named("mu_sigma") = h_ms,
                        Named("sigma_sigma") = h_ss);
}

// The expected information, by summing the observed Hessian against the mass
// over the finite support {0, ..., n}. The support being finite makes this an
// EXACT expectation rather than a quadrature, and the result does not depend
// on the data, only on the parameters.
// [[Rcpp::export]]
List betabinom_expected_hessian_cpp(NumericVector y, NumericVector mu,
                                    NumericVector sigma, double size) {
    int n = y.size();
    NumericVector h_mm(n), h_ms(n), h_ss(n);
    bool mu_s = (mu.size() == 1), si_s = (sigma.size() == 1);
    int N = (int) size;

    double last_m = R_NegInf, last_s = R_NegInf;
    double e_mm = 0, e_ms = 0, e_ss = 0;

    for (int i = 0; i < n; i++) {
        double m = mu_s ? mu[0] : mu[i];
        double s = si_s ? sigma[0] : sigma[i];

        if (m != last_m || s != last_s) {
            BBmap mp = bb_map(m, s);
            double lb0 = R::lbeta(mp.A, mp.B);
            e_mm = 0; e_ms = 0; e_ss = 0;
            for (int k = 0; k <= N; k++) {
                double p = std::exp(R::lchoose(size, k)
                                    + R::lbeta(k + mp.A, size - k + mp.B) - lb0);
                BBderiv d = bb_shape_derivs(k, size, mp.A, mp.B);
                e_mm += p * (d.lAA * mp.Am * mp.Am
                             + 2.0 * d.lAB * mp.Am * mp.Bm
                             + d.lBB * mp.Bm * mp.Bm
                             + d.lA * mp.Amm + d.lB * mp.Bmm);
                e_ms += p * (d.lAA * mp.Am * mp.As
                             + d.lAB * (mp.Am * mp.Bs + mp.As * mp.Bm)
                             + d.lBB * mp.Bm * mp.Bs
                             + d.lA * mp.Ams + d.lB * mp.Bms);
                e_ss += p * (d.lAA * mp.As * mp.As
                             + 2.0 * d.lAB * mp.As * mp.Bs
                             + d.lBB * mp.Bs * mp.Bs
                             + d.lA * mp.Ass + d.lB * mp.Bss);
            }
            last_m = m; last_s = s;
        }
        h_mm[i] = e_mm; h_ms[i] = e_ms; h_ss[i] = e_ss;
    }
    return List::create(Named("mu_mu") = h_mm,
                        Named("mu_sigma") = h_ms,
                        Named("sigma_sigma") = h_ss);
}

// [[Rcpp::export]]
NumericVector betabinom_logpmf_cpp(NumericVector y, NumericVector mu,
                                   NumericVector sigma, double size) {
    int n = y.size();
    NumericVector out(n);
    bool mu_s = (mu.size() == 1), si_s = (sigma.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = mu_s ? mu[0] : mu[i];
        double s = si_s ? sigma[0] : sigma[i];
        double A = m / s, B = (1.0 - m) / s;
        if (y[i] < 0 || y[i] > size || y[i] != std::floor(y[i])) {
            out[i] = R_NegInf;
        } else {
            out[i] = R::lchoose(size, y[i])
                   + R::lbeta(y[i] + A, size - y[i] + B) - R::lbeta(A, B);
        }
    }
    return out;
}
