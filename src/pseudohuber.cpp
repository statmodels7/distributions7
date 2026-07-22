#include <Rcpp.h>
using namespace Rcpp;

// Ratios of Bessel K derivatives to K_1, computed with exponentially scaled
// Bessel functions so the e^{-x} factors cancel and large arguments do not
// underflow.
//   r1 = K_1'(x) / K_1(x) = -(K_0 + K_2) / (2 K_1)
//   r2 = K_1''(x) / K_1(x) = (3 K_1 + K_3) / (4 K_1)
static void bessel_ratios(double x, double &r1, double &r2) {
    double k0 = R::bessel_k(x, 0.0, 2.0);
    double k1 = R::bessel_k(x, 1.0, 2.0);
    double k2 = R::bessel_k(x, 2.0, 2.0);
    double k3 = R::bessel_k(x, 3.0, 2.0);
    r1 = -(k0 + k2) / (2.0 * k1);
    r2 = (3.0 * k1 + k3) / (4.0 * k1);
}

// [[Rcpp::export]]
List pseudohuber_gradient_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu) {
    int n = y.size();
    NumericVector grad_mu(n), grad_sigma(n), grad_nu(n);

    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);
    bool nu_is_scalar = (nu.size() == 1);
    bool all_scalar = mu_is_scalar && sigma_is_scalar && nu_is_scalar;

    double m = 0, s = 0, v = 0, s2 = 0, sq_v = 0, r1 = 0, r2 = 0;

    if (all_scalar) {
        m = mu[0]; s = sigma[0]; v = nu[0];
        s2 = s * s;
        sq_v = std::sqrt(v);
        bessel_ratios(sq_v, r1, r2);
    }

    for(int i = 0; i < n; i++) {
        if (!all_scalar) {
            m = mu_is_scalar ? mu[0] : mu[i];
            s = sigma_is_scalar ? sigma[0] : sigma[i];
            v = nu_is_scalar ? nu[0] : nu[i];
            s2 = s * s;
            sq_v = std::sqrt(v);
            bessel_ratios(sq_v, r1, r2);
        }

        double res = y[i] - m;
        double res2 = res * res;
        double D = std::sqrt(v + res2 / s2);

        grad_mu[i] = res / (s2 * D);
        grad_sigma[i] = (res2 / (s2 * D) - 1.0) / s;
        grad_nu[i] = -0.5 * (1.0 / v + 1.0 / D + r1 / sq_v);
    }

    return List::create(Named("mu") = grad_mu, Named("sigma") = grad_sigma, Named("nu") = grad_nu);
}

// [[Rcpp::export]]
List pseudohuber_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu) {
    int n = y.size();
    NumericVector hess_mu_mu(n), hess_sigma_sigma(n), hess_nu_nu(n);
    NumericVector hess_mu_sigma(n), hess_mu_nu(n), hess_sigma_nu(n);

    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);
    bool nu_is_scalar = (nu.size() == 1);
    bool all_scalar = mu_is_scalar && sigma_is_scalar && nu_is_scalar;

    double m = 0, s = 0, v = 0, s2 = 0, s4 = 0, sq_v = 0, r1 = 0, r2 = 0, nu_const = 0;

    if (all_scalar) {
        m = mu[0]; s = sigma[0]; v = nu[0];
        s2 = s * s; s4 = s2 * s2;
        sq_v = std::sqrt(v);
        bessel_ratios(sq_v, r1, r2);
        nu_const = 0.5 / (v * v) + 0.25 * (r1 / (v * sq_v) + r1 * r1 / v - r2 / v);
    }

    for(int i = 0; i < n; i++) {
        if (!all_scalar) {
            m = mu_is_scalar ? mu[0] : mu[i];
            s = sigma_is_scalar ? sigma[0] : sigma[i];
            v = nu_is_scalar ? nu[0] : nu[i];
            s2 = s * s; s4 = s2 * s2;
            sq_v = std::sqrt(v);
            bessel_ratios(sq_v, r1, r2);
            nu_const = 0.5 / (v * v) + 0.25 * (r1 / (v * sq_v) + r1 * r1 / v - r2 / v);
        }

        double res = y[i] - m;
        double res2 = res * res;
        double D = std::sqrt(v + res2 / s2);
        double D3 = D * D * D;

        hess_mu_mu[i] = -v / (s2 * D3);
        hess_sigma_sigma[i] = (s4 - 3.0 * s2 * res2 / D + res2 * res2 / D3) / (s4 * s2);
        hess_nu_nu[i] = 0.25 / D3 + nu_const;

        hess_mu_sigma[i] = (-2.0 * v * s2 * res - res2 * res) / (s2 * std::pow(v * s2 + res2, 1.5));
        hess_mu_nu[i] = -res / (2.0 * s2 * D3);
        hess_sigma_nu[i] = -res2 / (2.0 * s2 * s * D3);
    }

    return List::create(
        Named("mu_mu") = hess_mu_mu, Named("sigma_sigma") = hess_sigma_sigma, Named("nu_nu") = hess_nu_nu,
        Named("mu_sigma") = hess_mu_sigma, Named("mu_nu") = hess_mu_nu, Named("sigma_nu") = hess_sigma_nu
    );
}
