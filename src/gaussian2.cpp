#include <Rcpp.h>
using namespace Rcpp;

// Gaussian in the mean and the VARIANCE, l = -log(2 pi v)/2 - r^2/(2 v),
// with r = y - mu. Every derivative below is written out from that; the
// expected ones use E[r] = 0, E[r^2] = v, E[r^4] = 3 v^2.

// [[Rcpp::export]]
List gaussian2_gradient_cpp(NumericVector y, NumericVector mu, NumericVector sigma2) {
    int n = y.size();
    NumericVector g_mu(n), g_v(n);
    bool m_s = (mu.size() == 1), v_s = (sigma2.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = m_s ? mu[0] : mu[i];
        double v = v_s ? sigma2[0] : sigma2[i];
        double r = y[i] - m;
        g_mu[i] = r / v;
        g_v[i] = (r * r / v - 1.0) / (2.0 * v);
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
        double r = y[i] - m, v2 = v * v, v3 = v2 * v;
        h_mm[i] = -1.0 / v;
        h_mv[i] = -r / v2;
        h_vv[i] = 1.0 / (2.0 * v2) - r * r / v3;
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
        h_mm[i] = -1.0 / v;
        h_mv[i] = 0.0;
        h_vv[i] = -1.0 / (2.0 * v * v);
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
        double r = y[i] - m, v2 = v * v, v3 = v2 * v, v4 = v2 * v2;
        a[i] = 0.0;
        b[i] = 1.0 / v2;
        c[i] = 2.0 * r / v3;
        d[i] = -1.0 / v3 + 3.0 * r * r / v4;
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
        double v2 = v * v, v3 = v2 * v;
        a[i] = 0.0;
        b[i] = 1.0 / v2;
        c[i] = 0.0;
        d[i] = 2.0 / v3;
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
        double r = y[i] - m, v2 = v * v, v3 = v2 * v, v4 = v2 * v2, v5 = v4 * v;
        a[i] = 0.0;
        b[i] = 0.0;
        c[i] = -2.0 / v3;
        d[i] = -6.0 * r / v4;
        e[i] = 3.0 / v4 - 12.0 * r * r / v5;
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
        double v2 = v * v, v3 = v2 * v, v4 = v2 * v2;
        a[i] = 0.0;
        b[i] = 0.0;
        c[i] = -2.0 / v3;
        d[i] = 0.0;
        e[i] = -9.0 / v4;
    }
    return List::create(Named("mu_mu_mu_mu") = a, Named("mu_mu_mu_sigma2") = b,
                        Named("mu_mu_sigma2_sigma2") = c,
                        Named("mu_sigma2_sigma2_sigma2") = d,
                        Named("sigma2_sigma2_sigma2_sigma2") = e);
}
