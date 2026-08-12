#include <Rcpp.h>
using namespace Rcpp;

// Gaussian in the mean and the VARIANCE, l = -log(2 pi v)/2 - r^2/(2 v),
// with r = y - mu. Every derivative below is written out from that; the
// expected ones use E[r] = 0, E[r^2] = v, E[r^4] = 3 v^2.
//
// The components are carried in u = 1/v and z2 = r^2/v rather than in powers
// of the variance, for the reason given in gaussian.cpp: a denominator
// carrying v^5 overflows at v = 1e61, and a component computed as a
// difference of two such ratios then loses the term that dominates it.

// [[Rcpp::export]]
List gaussian2_gradient_cpp(NumericVector y, NumericVector mu, NumericVector sigma2) {
    int n = y.size();
    NumericVector g_mu(n), g_v(n);
    bool m_s = (mu.size() == 1), v_s = (sigma2.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = m_s ? mu[0] : mu[i];
        double v = v_s ? sigma2[0] : sigma2[i];
        double u = 1.0 / v;
        double r = y[i] - m;
        double z2 = r * r / v;
        g_mu[i] = r * u;
        g_v[i] = 0.5 * (z2 - 1.0) * u;
    }
    return List::create(Named("mu") = g_mu, Named("sigma2") = g_v);
}

// [[Rcpp::export]]
List gaussian2_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma2) {
    int n = y.size();
    NumericVector h_mm(n), h_mv(n), h_vv(n);
    bool m_s = (mu.size() == 1), v_s = (sigma2.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = m_s ? mu[0] : mu[i];
        double v = v_s ? sigma2[0] : sigma2[i];
        double u = 1.0 / v, u2 = u * u;
        double r = y[i] - m;
        double z2 = r * r / v;
        h_mm[i] = -u;
        h_mv[i] = -r * u2;
        h_vv[i] = (0.5 - z2) * u2;
    }
    return List::create(Named("mu_mu") = h_mm, Named("mu_sigma2") = h_mv,
                        Named("sigma2_sigma2") = h_vv);
}

// [[Rcpp::export]]
List gaussian2_expected_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma2) {
    int n = y.size();
    NumericVector h_mm(n), h_mv(n), h_vv(n);
    bool v_s = (sigma2.size() == 1);

    for (int i = 0; i < n; i++) {
        double v = v_s ? sigma2[0] : sigma2[i];
        double u = 1.0 / v;
        h_mm[i] = -u;
        h_mv[i] = 0.0;
        h_vv[i] = -0.5 * u * u;
    }
    return List::create(Named("mu_mu") = h_mm, Named("mu_sigma2") = h_mv,
                        Named("sigma2_sigma2") = h_vv);
}

// [[Rcpp::export]]
List gaussian2_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector sigma2) {
    int n = y.size();
    NumericVector a(n), b(n), c(n), d(n);
    bool m_s = (mu.size() == 1), v_s = (sigma2.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = m_s ? mu[0] : mu[i];
        double v = v_s ? sigma2[0] : sigma2[i];
        double u = 1.0 / v, u2 = u * u, u3 = u2 * u;
        double r = y[i] - m;
        double z2 = r * r / v;
        a[i] = 0.0;
        b[i] = u2;
        c[i] = 2.0 * r * u3;
        d[i] = (3.0 * z2 - 1.0) * u3;
    }
    return List::create(Named("mu_mu_mu") = a, Named("mu_mu_sigma2") = b,
                        Named("mu_sigma2_sigma2") = c,
                        Named("sigma2_sigma2_sigma2") = d);
}

// [[Rcpp::export]]
List gaussian2_deriv3_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma2) {
    int n = y.size();
    NumericVector a(n), b(n), c(n), d(n);
    bool v_s = (sigma2.size() == 1);

    for (int i = 0; i < n; i++) {
        double v = v_s ? sigma2[0] : sigma2[i];
        double u = 1.0 / v, u2 = u * u;
        a[i] = 0.0;
        b[i] = u2;
        c[i] = 0.0;
        d[i] = 2.0 * u2 * u;
    }
    return List::create(Named("mu_mu_mu") = a, Named("mu_mu_sigma2") = b,
                        Named("mu_sigma2_sigma2") = c,
                        Named("sigma2_sigma2_sigma2") = d);
}

// [[Rcpp::export]]
List gaussian2_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector sigma2) {
    int n = y.size();
    NumericVector a(n), b(n), c(n), d(n), e(n);
    bool m_s = (mu.size() == 1), v_s = (sigma2.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = m_s ? mu[0] : mu[i];
        double v = v_s ? sigma2[0] : sigma2[i];
        double u = 1.0 / v, u2 = u * u, u3 = u2 * u, u4 = u2 * u2;
        double r = y[i] - m;
        double z2 = r * r / v;
        a[i] = 0.0;
        b[i] = 0.0;
        c[i] = -2.0 * u3;
        d[i] = -6.0 * r * u4;
        e[i] = (3.0 - 12.0 * z2) * u4;
    }
    return List::create(Named("mu_mu_mu_mu") = a, Named("mu_mu_mu_sigma2") = b,
                        Named("mu_mu_sigma2_sigma2") = c,
                        Named("mu_sigma2_sigma2_sigma2") = d,
                        Named("sigma2_sigma2_sigma2_sigma2") = e);
}

// [[Rcpp::export]]
List gaussian2_deriv4_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma2) {
    int n = y.size();
    NumericVector a(n), b(n), c(n), d(n), e(n);
    bool v_s = (sigma2.size() == 1);

    for (int i = 0; i < n; i++) {
        double v = v_s ? sigma2[0] : sigma2[i];
        double u = 1.0 / v, u2 = u * u, u3 = u2 * u, u4 = u2 * u2;
        a[i] = 0.0;
        b[i] = 0.0;
        c[i] = -2.0 * u3;
        d[i] = 0.0;
        e[i] = -9.0 * u4;
    }
    return List::create(Named("mu_mu_mu_mu") = a, Named("mu_mu_mu_sigma2") = b,
                        Named("mu_mu_sigma2_sigma2") = c,
                        Named("mu_sigma2_sigma2_sigma2") = d,
                        Named("sigma2_sigma2_sigma2_sigma2") = e);
}
