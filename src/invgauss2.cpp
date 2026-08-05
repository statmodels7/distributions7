#include <Rcpp.h>
using namespace Rcpp;

// Inverse gaussian in the mean and the SHAPE lambda, the classical
// parametrisation, with Var = mu^3 / lambda:
//
//   l = log(lambda)/2 - log(2 pi y^3)/2 - lambda (y-mu)^2 / (2 mu^2 y).
//
// The quadratic form collapses to -lambda (y-mu)^2/(2 mu^2 y) whose mu
// derivative is lambda (y-mu)/mu^3, so every order is a rational function of
// mu with y appearing linearly. The expected versions therefore need only
// E[y] = mu.

// [[Rcpp::export]]
List invgauss2_gradient_cpp(NumericVector y, NumericVector mu, NumericVector lambda) {
    int n = y.size();
    NumericVector g_mu(n), g_l(n);
    bool m_s = (mu.size() == 1), l_s = (lambda.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = m_s ? mu[0] : mu[i];
        double L = l_s ? lambda[0] : lambda[i];
        double r = y[i] - m, m2 = m * m, m3 = m2 * m;
        g_mu[i] = L * r / m3;
        g_l[i] = 0.5 / L - r * r / (2.0 * m2 * y[i]);
    }
    return List::create(Named("mu") = g_mu, Named("lambda") = g_l);
}

// [[Rcpp::export]]
List invgauss2_hessian_cpp(NumericVector y, NumericVector mu, NumericVector lambda) {
    int n = y.size();
    NumericVector h_mm(n), h_ml(n), h_ll(n);
    bool m_s = (mu.size() == 1), l_s = (lambda.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = m_s ? mu[0] : mu[i];
        double L = l_s ? lambda[0] : lambda[i];
        double m2 = m * m, m3 = m2 * m, m4 = m2 * m2;
        h_mm[i] = L * (2.0 * m - 3.0 * y[i]) / m4;
        h_ml[i] = (y[i] - m) / m3;
        h_ll[i] = -0.5 / (L * L);
    }
    return List::create(Named("mu_mu") = h_mm, Named("mu_lambda") = h_ml,
                        Named("lambda_lambda") = h_ll);
}

// [[Rcpp::export]]
List invgauss2_expected_hessian_cpp(NumericVector y, NumericVector mu, NumericVector lambda) {
    int n = y.size();
    NumericVector h_mm(n), h_ml(n), h_ll(n);
    bool m_s = (mu.size() == 1), l_s = (lambda.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = m_s ? mu[0] : mu[i];
        double L = l_s ? lambda[0] : lambda[i];
        h_mm[i] = -L / (m * m * m);
        h_ml[i] = 0.0;
        h_ll[i] = -0.5 / (L * L);
    }
    return List::create(Named("mu_mu") = h_mm, Named("mu_lambda") = h_ml,
                        Named("lambda_lambda") = h_ll);
}

// [[Rcpp::export]]
List invgauss2_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector lambda,
                          bool expected) {
    int n = y.size();
    NumericVector a(n), b(n), c(n), d(n);
    bool m_s = (mu.size() == 1), l_s = (lambda.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = m_s ? mu[0] : mu[i];
        double L = l_s ? lambda[0] : lambda[i];
        double yy = expected ? m : y[i];
        double m2 = m * m, m4 = m2 * m2, m5 = m4 * m;
        a[i] = L * (12.0 * yy - 6.0 * m) / m5;
        b[i] = (2.0 * m - 3.0 * yy) / m4;
        c[i] = 0.0;
        d[i] = 1.0 / (L * L * L);
    }
    return List::create(Named("mu_mu_mu") = a, Named("mu_mu_lambda") = b,
                        Named("mu_lambda_lambda") = c,
                        Named("lambda_lambda_lambda") = d);
}

// [[Rcpp::export]]
List invgauss2_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector lambda,
                          bool expected) {
    int n = y.size();
    NumericVector a(n), b(n), c(n), d(n), e(n);
    bool m_s = (mu.size() == 1), l_s = (lambda.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = m_s ? mu[0] : mu[i];
        double L = l_s ? lambda[0] : lambda[i];
        double yy = expected ? m : y[i];
        double m2 = m * m, m4 = m2 * m2, m5 = m4 * m, m6 = m5 * m;
        double L2 = L * L;
        a[i] = L * (24.0 * m - 60.0 * yy) / m6;
        b[i] = (12.0 * yy - 6.0 * m) / m5;
        c[i] = 0.0;
        d[i] = 0.0;
        e[i] = -3.0 / (L2 * L2);
    }
    return List::create(Named("mu_mu_mu_mu") = a, Named("mu_mu_mu_lambda") = b,
                        Named("mu_mu_lambda_lambda") = c,
                        Named("mu_lambda_lambda_lambda") = d,
                        Named("lambda_lambda_lambda_lambda") = e);
}
