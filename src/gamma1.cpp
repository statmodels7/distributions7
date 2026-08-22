#include <Rcpp.h>
#include "d7_par.h"
#include "psi_diff.h"
using namespace Rcpp;

// Gamma in the mean and the DISPERSION, Var = phi mu^2, which is the GLM
// parametrization. Writing s = 1/phi for the shape and z = y/mu,
//
//   l = s log s - s log mu - lgamma(s) + (s-1) log y - s z.
//
// Every derivative in phi goes through s, so they are the derivatives in s
// carried across by the one-variable chain rule with
//   s' = -s^2, s'' = 2 s^3, s''' = -6 s^4, s'''' = 24 s^5.
// The expected versions use E[z] = 1 and E[log z] = digamma(s) - log s, the
// second of which is what makes E[l_s] vanish.
//
// Each kernel is a per-observation map -- observation i's derivatives are
// written to slot i and nothing else -- so the loop runs through
// d7::par_for, whose decomposition over output elements keeps the result
// bit-identical at any thread count.

namespace {

struct Fs {
    double f1, f2, f3, f4;      // derivatives in s
    double m1, m2, m3, m4;      // derivatives in mu
    double ms1, mms1, mmms1;    // mixed, one s
};

inline Fs gamma1_parts(double y, double m, double phi) {
    double s = 1.0 / phi;
    double z = y / m;
    double m2 = m * m, m3 = m2 * m, m4 = m2 * m2;
    Fs o;
    // [log s - psi(s)] + [log z - (z-1)]: two cancelling pairs, each
    // written out.  See psi_diff.h.
    o.f1 = d7::psi_log_rest(s) + d7::psi_Ew(z - 1.0);
    o.f2 = d7::psi1_rest(s);
    o.f3 = d7::psi2_rest(s);
    o.f4 = d7::psi3_rest(s);
    o.m1 = s * (z - 1.0) / m;
    o.m2 = s * (1.0 - 2.0 * z) / m2;
    o.m3 = s * (6.0 * z - 2.0) / m3;
    o.m4 = s * (6.0 - 24.0 * z) / m4;
    o.ms1 = (z - 1.0) / m;
    o.mms1 = (1.0 - 2.0 * z) / m2;
    o.mmms1 = (6.0 * z - 2.0) / m3;
    return o;
}

// The expected versions of the same quantities, at E[z] = 1 and
// E[log z] = digamma(s) - log s.
inline Fs gamma1_parts_expected(double m, double phi) {
    double s = 1.0 / phi;
    double m2 = m * m, m3 = m2 * m, m4 = m2 * m2;
    Fs o;
    o.f1 = 0.0;
    o.f2 = d7::psi1_rest(s);
    o.f3 = d7::psi2_rest(s);
    o.f4 = d7::psi3_rest(s);
    o.m1 = 0.0;
    o.m2 = -s / m2;
    o.m3 = 4.0 * s / m3;
    o.m4 = -18.0 * s / m4;
    o.ms1 = 0.0;
    o.mms1 = -1.0 / m2;
    o.mmms1 = 4.0 / m3;
    return o;
}

} // namespace

// [[Rcpp::export]]
List gamma1_gradient_cpp(NumericVector y, NumericVector mu, NumericVector phi,
                         int threads = 1) {
    int n = y.size();
    NumericVector g_mu(n), g_ph(n);
    bool m_s = (mu.size() == 1), p_s = (phi.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *pp = phi.begin();
    double *gm = g_mu.begin(), *gp = g_ph.begin();

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m_s ? mp[0] : mp[i];
        double p = p_s ? pp[0] : pp[i];
        double s = 1.0 / p;
        Fs o = gamma1_parts(yp[i], m, p);
        gm[i] = o.m1;
        gp[i] = o.f1 * (-s * s);
    });
    return List::create(Named("mu") = g_mu, Named("phi") = g_ph);
}

// [[Rcpp::export]]
List gamma1_hessian_cpp(NumericVector y, NumericVector mu, NumericVector phi,
                        int threads = 1) {
    int n = y.size();
    NumericVector h_mm(n), h_mp(n), h_pp(n);
    bool m_s = (mu.size() == 1), p_s = (phi.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *pp = phi.begin();
    double *hmm = h_mm.begin(), *hmp = h_mp.begin(), *hpp = h_pp.begin();

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m_s ? mp[0] : mp[i];
        double p = p_s ? pp[0] : pp[i];
        double s = 1.0 / p, s1 = -s * s, s2 = 2.0 * s * s * s;
        Fs o = gamma1_parts(yp[i], m, p);
        hmm[i] = o.m2;
        hmp[i] = o.ms1 * s1;
        hpp[i] = o.f2 * s1 * s1 + o.f1 * s2;
    });
    return List::create(Named("mu_mu") = h_mm, Named("mu_phi") = h_mp,
                        Named("phi_phi") = h_pp);
}

// [[Rcpp::export]]
List gamma1_expected_hessian_cpp(NumericVector y, NumericVector mu,
                                 NumericVector phi, int threads = 1) {
    int n = y.size();
    NumericVector h_mm(n), h_mp(n), h_pp(n);
    bool m_s = (mu.size() == 1), p_s = (phi.size() == 1);
    const double *mp = mu.begin(), *pp = phi.begin();
    double *hmm = h_mm.begin(), *hmp = h_mp.begin(), *hpp = h_pp.begin();

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m_s ? mp[0] : mp[i];
        double p = p_s ? pp[0] : pp[i];
        double s = 1.0 / p, s1 = -s * s;
        Fs o = gamma1_parts_expected(m, p);
        hmm[i] = o.m2;
        hmp[i] = 0.0;
        hpp[i] = o.f2 * s1 * s1;
    });
    return List::create(Named("mu_mu") = h_mm, Named("mu_phi") = h_mp,
                        Named("phi_phi") = h_pp);
}

namespace {

inline void gamma1_d3(const Fs &o, double s, double *out) {
    double s1 = -s * s, s2 = 2.0 * s * s * s, s3 = -6.0 * s * s * s * s;
    out[0] = o.m3;                                   // mu mu mu
    out[1] = o.mms1 * s1;                            // mu mu phi
    out[2] = o.ms1 * s2;                             // mu phi phi
    out[3] = o.f3 * s1 * s1 * s1 + 3.0 * o.f2 * s1 * s2 + o.f1 * s3;
}

inline void gamma1_d4(const Fs &o, double s, double *out) {
    double s1 = -s * s, s2 = 2.0 * s * s * s, s3 = -6.0 * s * s * s * s,
           s4 = 24.0 * s * s * s * s * s;
    out[0] = o.m4;                                   // mu mu mu mu
    out[1] = o.mmms1 * s1;                           // mu mu mu phi
    out[2] = o.mms1 * s2;                            // mu mu phi phi
    out[3] = o.ms1 * s3;                             // mu phi phi phi
    out[4] = o.f4 * s1 * s1 * s1 * s1
           + 6.0 * o.f3 * s1 * s1 * s2
           + 3.0 * o.f2 * s2 * s2
           + 4.0 * o.f2 * s1 * s3
           + o.f1 * s4;
}

} // namespace

// [[Rcpp::export]]
List gamma1_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector phi,
                       int threads = 1) {
    int n = y.size();
    NumericVector a(n), b(n), c(n), d(n);
    bool m_s = (mu.size() == 1), p_s = (phi.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *pp = phi.begin();
    double *ap = a.begin(), *bp = b.begin(), *cp = c.begin(), *dp = d.begin();

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m_s ? mp[0] : mp[i];
        double p = p_s ? pp[0] : pp[i];
        double v[4];
        gamma1_d3(gamma1_parts(yp[i], m, p), 1.0 / p, v);
        ap[i] = v[0]; bp[i] = v[1]; cp[i] = v[2]; dp[i] = v[3];
    });
    return List::create(Named("mu_mu_mu") = a, Named("mu_mu_phi") = b,
                        Named("mu_phi_phi") = c, Named("phi_phi_phi") = d);
}

// [[Rcpp::export]]
List gamma1_deriv3_expected_cpp(NumericVector y, NumericVector mu,
                                NumericVector phi, int threads = 1) {
    int n = y.size();
    NumericVector a(n), b(n), c(n), d(n);
    bool m_s = (mu.size() == 1), p_s = (phi.size() == 1);
    const double *mp = mu.begin(), *pp = phi.begin();
    double *ap = a.begin(), *bp = b.begin(), *cp = c.begin(), *dp = d.begin();

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m_s ? mp[0] : mp[i];
        double p = p_s ? pp[0] : pp[i];
        double v[4];
        gamma1_d3(gamma1_parts_expected(m, p), 1.0 / p, v);
        ap[i] = v[0]; bp[i] = v[1]; cp[i] = v[2]; dp[i] = v[3];
    });
    return List::create(Named("mu_mu_mu") = a, Named("mu_mu_phi") = b,
                        Named("mu_phi_phi") = c, Named("phi_phi_phi") = d);
}

// [[Rcpp::export]]
List gamma1_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector phi,
                       int threads = 1) {
    int n = y.size();
    NumericVector a(n), b(n), c(n), d(n), e(n);
    bool m_s = (mu.size() == 1), p_s = (phi.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *pp = phi.begin();
    double *ap = a.begin(), *bp = b.begin(), *cp = c.begin(), *dp = d.begin(),
           *ep = e.begin();

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m_s ? mp[0] : mp[i];
        double p = p_s ? pp[0] : pp[i];
        double v[5];
        gamma1_d4(gamma1_parts(yp[i], m, p), 1.0 / p, v);
        ap[i] = v[0]; bp[i] = v[1]; cp[i] = v[2]; dp[i] = v[3]; ep[i] = v[4];
    });
    return List::create(Named("mu_mu_mu_mu") = a, Named("mu_mu_mu_phi") = b,
                        Named("mu_mu_phi_phi") = c, Named("mu_phi_phi_phi") = d,
                        Named("phi_phi_phi_phi") = e);
}

// [[Rcpp::export]]
List gamma1_deriv4_expected_cpp(NumericVector y, NumericVector mu,
                                NumericVector phi, int threads = 1) {
    int n = y.size();
    NumericVector a(n), b(n), c(n), d(n), e(n);
    bool m_s = (mu.size() == 1), p_s = (phi.size() == 1);
    const double *mp = mu.begin(), *pp = phi.begin();
    double *ap = a.begin(), *bp = b.begin(), *cp = c.begin(), *dp = d.begin(),
           *ep = e.begin();

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m_s ? mp[0] : mp[i];
        double p = p_s ? pp[0] : pp[i];
        double v[5];
        gamma1_d4(gamma1_parts_expected(m, p), 1.0 / p, v);
        ap[i] = v[0]; bp[i] = v[1]; cp[i] = v[2]; dp[i] = v[3]; ep[i] = v[4];
    });
    return List::create(Named("mu_mu_mu_mu") = a, Named("mu_mu_mu_phi") = b,
                        Named("mu_mu_phi_phi") = c, Named("mu_phi_phi_phi") = d,
                        Named("phi_phi_phi_phi") = e);
}
