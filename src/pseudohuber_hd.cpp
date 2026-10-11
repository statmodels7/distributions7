#include <Rcpp.h>
#include <cmath>
#include "pt_pseudohuber.h"
#include "pt_sqrt.h"
using namespace Rcpp;

// Observed third/fourth-order derivatives of the Pseudo-Huber log-density,
// transcribed from the Wolfram output. Everything is rewritten in terms of
//   r = mu - y,   S = nu*sigma^2 + r^2
// so that (nu + r^2/sigma^2)^(k/2) = S^(k/2) / sigma^k.
//
// Bessel functions appear only in the pure-nu derivatives, from numericals7's
// n7_bessel_k (pt_bessel_k.h), which returns R's values (the normalizing
// constant depends on sigma and nu separably). Every such term is homogeneous of
// the same degree in K over the same degree in K1, so the exponentially scaled
// Bessel functions may be used throughout: the e^{-x} factors cancel exactly and
// large nu no longer overflows.
//
// The expected higher-order derivatives have no closed form (the Wolfram
// integrals do not converge symbolically); they are obtained through the
// `approx` machinery in R.

// [[Rcpp::export]]
List pseudohuber_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_sigma(n), mu_mu_nu(n), mu_sigma_sigma(n), mu_sigma_nu(n),
                  mu_nu_nu(n), sigma_sigma_sigma(n), sigma_sigma_nu(n), sigma_nu_nu(n), nu_nu_nu(n);
    bool mu_s = (mu.size() == 1), sig_s = (sigma.size() == 1), nu_s = (nu.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = mu_s ? mu[0] : mu[i];
        double s = sig_s ? sigma[0] : sigma[i];
        double v = nu_s ? nu[0] : nu[i];

        double s2 = s * s, s3 = s2 * s;
        double r = m - y[i];
        const d7::PhA PA = d7::ph_a(r, s, v);
        const double A = PA.A, p = PA.rho, p2 = p * p, k = PA.kap2, q = s / A;
        double sv = d7::sqrt_cr(v);

        double k0 = d7::bessel_k_scaled(sv, 0.0);
        double k1 = d7::bessel_k_scaled(sv, 1.0);
        double k2 = d7::bessel_k_scaled(sv, 2.0);
        double k3 = d7::bessel_k_scaled(sv, 3.0);
        double k4 = d7::bessel_k_scaled(sv, 4.0);

        mu_mu_mu[i] = d7::pseudohuber_d3_mu_mu_mu(s, PA);
        mu_mu_sigma[i] = k * (2.0 * k - p2) / (s2 * A);
        mu_mu_nu[i] = q * (k - 2.0 * p2) / (2.0 * A * A);
        mu_sigma_sigma[i] = -p * (6.0 - 7.0 * k * p2 - 4.0 * p2 * p2) / s3;
        mu_sigma_nu[i] = p * (p2 - 2.0 * k) / (2.0 * A * A);
        mu_nu_nu[i] = -3.0 * q * q * q * p / (4.0 * A);
        sigma_sigma_sigma[i] = d7::pseudohuber_d3_sigma_sigma_sigma(r, s, PA);
        sigma_sigma_nu[i] = 3.0 * k * p2 / (2.0 * s * A);
        sigma_nu_nu[i] = 3.0 * q * q * p2 / (4.0 * A);
        const d7::PhNu P = d7::ph_nu(v);
        nu_nu_nu[i] = d7::pseudohuber_d3_nu_nu_nu(v, s, PA, P, k4);
    }

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu, Named("mu_mu_sigma") = mu_mu_sigma, Named("mu_mu_nu") = mu_mu_nu,
        Named("mu_sigma_sigma") = mu_sigma_sigma, Named("mu_sigma_nu") = mu_sigma_nu,
        Named("mu_nu_nu") = mu_nu_nu, Named("sigma_sigma_sigma") = sigma_sigma_sigma,
        Named("sigma_sigma_nu") = sigma_sigma_nu, Named("sigma_nu_nu") = sigma_nu_nu,
        Named("nu_nu_nu") = nu_nu_nu
    );
}

// [[Rcpp::export]]
List pseudohuber_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_sigma(n), mu_mu_mu_nu(n), mu_mu_sigma_sigma(n),
                  mu_mu_sigma_nu(n), mu_mu_nu_nu(n), mu_sigma_sigma_sigma(n), mu_sigma_sigma_nu(n),
                  mu_sigma_nu_nu(n), mu_nu_nu_nu(n), sigma_sigma_sigma_sigma(n),
                  sigma_sigma_sigma_nu(n), sigma_sigma_nu_nu(n), sigma_nu_nu_nu(n), nu_nu_nu_nu(n);
    bool mu_s = (mu.size() == 1), sig_s = (sigma.size() == 1), nu_s = (nu.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = mu_s ? mu[0] : mu[i];
        double s = sig_s ? sigma[0] : sigma[i];
        double v = nu_s ? nu[0] : nu[i];

        double s2 = s * s, s3 = s2 * s, s4 = s2 * s2, s5 = s4 * s;
        double r = m - y[i];
        const d7::PhA PA = d7::ph_a(r, s, v);
        const double A = PA.A, p = PA.rho, p2 = p * p, p4 = p2 * p2,
                     k = PA.kap2, q = s / A, A2 = A * A;
        double sv = d7::sqrt_cr(v), v2 = v * v, v3 = v2 * v, v4 = v2 * v2, v32 = v * sv;

        double k0 = d7::bessel_k_scaled(sv, 0.0);
        double k1 = d7::bessel_k_scaled(sv, 1.0);
        double k2 = d7::bessel_k_scaled(sv, 2.0);
        double k3 = d7::bessel_k_scaled(sv, 3.0);
        double k4 = d7::bessel_k_scaled(sv, 4.0);
        double k5 = d7::bessel_k_scaled(sv, 5.0);
        double k1_2 = k1 * k1, k1_3 = k1_2 * k1, k1_4 = k1_2 * k1_2;
        double k0_2 = k0 * k0, k0_3 = k0_2 * k0, k0_4 = k0_2 * k0_2;
        double k2_2 = k2 * k2, k2_3 = k2_2 * k2, k2_4 = k2_2 * k2_2;

        mu_mu_mu_mu[i] = 3.0 * k * (k - 4.0 * p2) / (s * A * A2);
        mu_mu_mu_sigma[i] = 3.0 * k * p * (p2 - 4.0 * k) / (s2 * A2);
        mu_mu_mu_nu[i] = 3.0 * q * p * (2.0 * p2 - 3.0 * k) / (2.0 * A * A2);
        mu_mu_sigma_sigma[i] = 3.0 * k * k * (3.0 * p2 - 2.0 * k) / (s3 * A);
        // numerator collapses to -2 r^4 + 11 nu sigma^2 r^2 - 2 nu^2 sigma^4
        mu_mu_sigma_nu[i] = (-2.0 * p4 + 11.0 * k * p2 - 2.0 * k * k) / (2.0 * A * A2);
        mu_mu_nu_nu[i] = -3.0 * q * q * q * (k - 4.0 * p2) / (4.0 * A2);
        mu_sigma_sigma_sigma[i] = -3.0 * p * (-8.0 + 16.0 * p2 - 15.0 * k * p4
            - 10.0 * p4 * p2) / s4;
        mu_sigma_sigma_nu[i] = 3.0 * k * p * (2.0 * k - 3.0 * p2) / (2.0 * s * A2);
        mu_sigma_nu_nu[i] = 3.0 * q * q * p * (2.0 * k - 3.0 * p2) / (4.0 * A2);
        mu_nu_nu_nu[i] = 15.0 * q * q * q * q * q * p / (8.0 * A);
        sigma_sigma_sigma_sigma[i] = 6.0 / s4 + (r / s5) * p *
            (-60.0 + 75.0 * p2 - 54.0 * p4 + 15.0 * p4 * p2);
        sigma_sigma_sigma_nu[i] = 3.0 * k * p2 * (p2 - 4.0 * k) / (2.0 * s2 * A);
        sigma_sigma_nu_nu[i] = 3.0 * q * p2 * (2.0 * p2 - 3.0 * k) / (4.0 * A2);
        sigma_nu_nu_nu[i] = -15.0 * q * q * q * q * p2 / (8.0 * A);

        // d^4 l/dnu^4 = (15/16) D^-7 - B''''(nu), with q = 1/D
        const d7::PhNu P = d7::ph_nu(v);
        const double q2 = q * q;
        nu_nu_nu_nu[i] = 0.9375 * q2 * q2 * q2 * q - P.dB[4];
    }

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu, Named("mu_mu_mu_sigma") = mu_mu_mu_sigma,
        Named("mu_mu_mu_nu") = mu_mu_mu_nu, Named("mu_mu_sigma_sigma") = mu_mu_sigma_sigma,
        Named("mu_mu_sigma_nu") = mu_mu_sigma_nu, Named("mu_mu_nu_nu") = mu_mu_nu_nu,
        Named("mu_sigma_sigma_sigma") = mu_sigma_sigma_sigma,
        Named("mu_sigma_sigma_nu") = mu_sigma_sigma_nu, Named("mu_sigma_nu_nu") = mu_sigma_nu_nu,
        Named("mu_nu_nu_nu") = mu_nu_nu_nu,
        Named("sigma_sigma_sigma_sigma") = sigma_sigma_sigma_sigma,
        Named("sigma_sigma_sigma_nu") = sigma_sigma_sigma_nu,
        Named("sigma_sigma_nu_nu") = sigma_sigma_nu_nu, Named("sigma_nu_nu_nu") = sigma_nu_nu_nu,
        Named("nu_nu_nu_nu") = nu_nu_nu_nu
    );
}
