#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Observed third/fourth-order derivatives of the Pseudo-Huber log-density,
// transcribed from the Wolfram output. Everything is rewritten in terms of
//   r = mu - y,   S = nu*sigma^2 + r^2
// so that (nu + r^2/sigma^2)^(k/2) = S^(k/2) / sigma^k.
//
// Bessel functions appear only in the pure-nu derivatives (the normalising
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

        double s2 = s * s, s3 = s2 * s, s4 = s2 * s2, s5 = s4 * s;
        double r = m - y[i], r2 = r * r, r4 = r2 * r2, r6 = r4 * r2;
        double S = v * s2 + r2;
        double S12 = std::sqrt(S), S32 = S * S12, S52 = S * S * S12;
        double sv = std::sqrt(v), v2 = v * v, v32 = v * sv;

        double k0 = R::bessel_k(sv, 0.0, 2.0);
        double k1 = R::bessel_k(sv, 1.0, 2.0);
        double k2 = R::bessel_k(sv, 2.0, 2.0);
        double k3 = R::bessel_k(sv, 3.0, 2.0);
        double k4 = R::bessel_k(sv, 4.0, 2.0);
        double A = k0 + k2;

        mu_mu_mu[i] = 3.0 * v * s * r / S52;
        mu_mu_sigma[i] = v * (2.0 * v * s2 - r2) / S52;
        mu_mu_nu[i] = s * (v * s2 - 2.0 * r2) / (2.0 * S52);
        mu_sigma_sigma[i] = (6.0 * S * S - 7.0 * v * s2 * r2 - 4.0 * r4) * (-r) / (s3 * S52);
        mu_sigma_nu[i] = (-2.0 * v * s2 + r2) * r / (2.0 * S52);
        mu_nu_nu[i] = -3.0 * s3 * r / (4.0 * S52);
        sigma_sigma_sigma[i] = -2.0 / s3 + 12.0 * r2 / (s4 * S12)
            - 9.0 * r4 / (s4 * S32) + 3.0 * r6 / (s4 * S52);
        sigma_sigma_nu[i] = 3.0 * v * s * r2 / (2.0 * S52);
        sigma_nu_nu[i] = 3.0 * s2 * r2 / (4.0 * S52);

        nu_nu_nu[i] = (
            -32.0 / (v2 * v)
            - 12.0 * s5 / S52
            + 6.0 * A / (v2 * sv * k1)
            - 3.0 * A * A / (v2 * k1 * k1)
            + A * A * A / (v32 * k1 * k1 * k1)
            - 3.0 * A * (3.0 * k1 + k3) / (2.0 * v32 * k1 * k1)
            + 2.0 * (3.0 + k3 / k1) / v2
            + (2.0 * (3.0 * k1 + k3) + sv * (3.0 * k0 + 4.0 * k2 + k4)) / (2.0 * v2 * k1)
        ) / 32.0;
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

        double s2 = s * s, s3 = s2 * s, s4 = s2 * s2, s5 = s4 * s, s7 = s5 * s2;
        double r = m - y[i], r2 = r * r, r4 = r2 * r2, r6 = r4 * r2, r8 = r4 * r4;
        double S = v * s2 + r2;
        double S12 = std::sqrt(S), S32 = S * S12, S52 = S * S * S12, S72 = S * S * S * S12;
        double sv = std::sqrt(v), v2 = v * v, v3 = v2 * v, v4 = v2 * v2, v32 = v * sv;

        double k0 = R::bessel_k(sv, 0.0, 2.0);
        double k1 = R::bessel_k(sv, 1.0, 2.0);
        double k2 = R::bessel_k(sv, 2.0, 2.0);
        double k3 = R::bessel_k(sv, 3.0, 2.0);
        double k4 = R::bessel_k(sv, 4.0, 2.0);
        double k5 = R::bessel_k(sv, 5.0, 2.0);
        double k1_2 = k1 * k1, k1_3 = k1_2 * k1, k1_4 = k1_2 * k1_2;
        double k0_2 = k0 * k0, k0_3 = k0_2 * k0, k0_4 = k0_2 * k0_2;
        double k2_2 = k2 * k2, k2_3 = k2_2 * k2, k2_4 = k2_2 * k2_2;

        mu_mu_mu_mu[i] = 3.0 * v * s * (v * s2 - 4.0 * r2) / S72;
        mu_mu_mu_sigma[i] = 3.0 * v * (-4.0 * v * s2 + r2) * r / S72;
        mu_mu_mu_nu[i] = 3.0 * s * (-3.0 * v * s2 + 2.0 * r2) * r / (2.0 * S72);
        mu_mu_sigma_sigma[i] = 3.0 * v2 * s * (-2.0 * v * s2 + 3.0 * r2) / S72;
        // numerator collapses to -2 r^4 + 11 nu sigma^2 r^2 - 2 nu^2 sigma^4
        mu_mu_sigma_nu[i] = (-2.0 * r4 + 11.0 * v * s2 * r2 - 2.0 * v2 * s4) / (2.0 * S72);
        mu_mu_nu_nu[i] = -3.0 * s3 * (v * s2 - 4.0 * r2) / (4.0 * S72);
        mu_sigma_sigma_sigma[i] = 3.0 * (-8.0 * S * S * S + 16.0 * S * S * r2
            - 15.0 * v * s2 * r4 - 10.0 * r6) * (-r) / (s4 * S72);
        mu_sigma_sigma_nu[i] = 3.0 * v * s * (2.0 * v * s2 - 3.0 * r2) * r / (2.0 * S72);
        mu_sigma_nu_nu[i] = 3.0 * s2 * (2.0 * v * s2 - 3.0 * r2) * r / (4.0 * S72);
        mu_nu_nu_nu[i] = 15.0 * s5 * r / (8.0 * S72);
        sigma_sigma_sigma_sigma[i] = 6.0 / s4 - 60.0 * r2 / (s5 * S12) + 75.0 * r4 / (s5 * S32)
            - 54.0 * r6 / (s5 * S52) + 15.0 * r8 / (s5 * S72);
        sigma_sigma_sigma_nu[i] = 3.0 * v * (-4.0 * v * s2 + r2) * r2 / (2.0 * S72);
        sigma_sigma_nu_nu[i] = 3.0 * s * (-3.0 * v * s2 + 2.0 * r2) * r2 / (4.0 * S72);
        sigma_nu_nu_nu[i] = -15.0 * s4 * r2 / (8.0 * S72);

        double N = 6.0 * v2 * k0_4
            + (768.0 + v * (-180.0 + 17.0 * v + 240.0 * v3 * s7 / S72)) * k1_4
            + 6.0 * v2 * k2_4
            + 24.0 * k0_3 * (-(v32 * k1) + v2 * k2)
            - 12.0 * k1 * (2.0 * v32 * k2_3 + v2 * k2_2 * k3)
            - 12.0 * v * k0_2 * ((-5.0 + 2.0 * v) * k1_2 - 3.0 * v * k2_2
                                 + k1 * (6.0 * sv * k2 + v * k3))
            - 4.0 * k0 * (6.0 * (5.0 - 3.0 * v) * sv * k1_3 - 6.0 * v2 * k2_3
                          + 6.0 * k1 * (3.0 * v32 * k2_2 + v2 * k2 * k3)
                          + v * k1_2 * ((-30.0 + 11.0 * v) * k2 - 9.0 * sv * k3 - v * k4))
            + v * k1_2 * (-20.0 * (-3.0 + v) * k2_2 + 3.0 * v * k3 * k3
                          + 4.0 * k2 * (9.0 * sv * k3 + v * k4))
            + k1_3 * (60.0 * (-2.0 + v) * sv * k2
                      + v * ((-60.0 + 13.0 * v) * k3 - 12.0 * sv * k4 - v * k5));

        nu_nu_nu_nu[i] = N / (256.0 * v4 * k1_4);
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
