#include <Rcpp.h>
#include "d7_par.h"
using namespace Rcpp;

// Third/fourth-order derivatives of the Gumbel (maxima) log-density. With
// z = (y - mu)/sigma and w = exp(-z) the log-likelihood is
// l = -log(sigma) - z - w, and the two generating rules
//   d/dmu    [P(z, w)/sigma^k] = (-P_z + w P_w) / sigma^(k+1)
//   d/dsigma [P(z, w)/sigma^k] = (-k P - z P_z + z w P_w) / sigma^(k+1)
// keep every derivative a polynomial in z and z^j w over a power of sigma.

// E[u (log u)^k] for a standard exponential u, i.e. Gamma^(k)(2); w = exp(-z)
// is standard exponential and z = -log(w), so E[z^k w] = (-1)^k Gamma^(k)(2)
// and E[z] = gamma (the Euler constant).
static void gumbel_gamma_moments(double g[5]) {
    const double euler = 0.577215664901532860606512090082;
    const double zeta3 = 1.202056903159594285399738161511;
    const double psi  = 1.0 - euler;
    const double psi1 = M_PI * M_PI / 6.0 - 1.0;
    const double psi2 = 2.0 - 2.0 * zeta3;
    const double psi3 = M_PI * M_PI * M_PI * M_PI / 15.0 - 6.0;
    g[0] = 1.0;
    g[1] = psi;
    g[2] = psi * psi + psi1;
    g[3] = psi * psi * psi + 3.0 * psi * psi1 + psi2;
    g[4] = psi * psi * psi * psi + 6.0 * psi * psi * psi1
         + 4.0 * psi * psi2 + 3.0 * psi1 * psi1 + psi3;
}

// [[Rcpp::export]]
List gumbel_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_sigma(n), mu_sigma_sigma(n), sigma_sigma_sigma(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double s3 = s * s * s;
        double z = (y[i] - m) / s;
        double w = std::exp(-z);
        double z2 = z * z;

        mu_mu_mu[i] = -w / s3;
        mu_mu_sigma[i] = w * (2.0 - z) / s3;
        mu_sigma_sigma[i] = (2.0 - 2.0 * w + 4.0 * z * w - z2 * w) / s3;
        sigma_sigma_sigma[i] = (-2.0 + 6.0 * z - 6.0 * z * w
                                + 6.0 * z2 * w - z2 * z * w) / s3;
    });

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_sigma") = mu_mu_sigma,
        Named("mu_sigma_sigma") = mu_sigma_sigma,
        Named("sigma_sigma_sigma") = sigma_sigma_sigma
    );
}

// [[Rcpp::export]]
List gumbel_deriv3_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_sigma(n), mu_sigma_sigma(n), sigma_sigma_sigma(n);
    bool sigma_is_scalar = (sigma.size() == 1);
    const double euler = 0.577215664901532860606512090082;
    double g[5];
    gumbel_gamma_moments(g);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double s3 = s * s * s;

        mu_mu_mu[i] = -1.0 / s3;
        mu_mu_sigma[i] = (3.0 - euler) / s3;
        mu_sigma_sigma[i] = (4.0 * euler - 4.0 - g[2]) / s3;
        sigma_sigma_sigma[i] = (4.0 + 6.0 * g[2] + g[3]) / s3;
    });

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_sigma") = mu_mu_sigma,
        Named("mu_sigma_sigma") = mu_sigma_sigma,
        Named("sigma_sigma_sigma") = sigma_sigma_sigma
    );
}

// [[Rcpp::export]]
List gumbel_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_sigma(n), mu_mu_sigma_sigma(n),
                  mu_sigma_sigma_sigma(n), sigma_sigma_sigma_sigma(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double s4 = s * s * s * s;
        double z = (y[i] - m) / s;
        double w = std::exp(-z);
        double z2 = z * z, z3 = z2 * z, z4 = z2 * z2;

        mu_mu_mu_mu[i] = -w / s4;
        mu_mu_mu_sigma[i] = w * (3.0 - z) / s4;
        mu_mu_sigma_sigma[i] = -w * (6.0 - 6.0 * z + z2) / s4;
        mu_sigma_sigma_sigma[i] = (-6.0 + 6.0 * w - 18.0 * z * w
                                   + 9.0 * z2 * w - z3 * w) / s4;
        sigma_sigma_sigma_sigma[i] = (6.0 - 24.0 * z + 24.0 * z * w
                                      - 36.0 * z2 * w + 12.0 * z3 * w
                                      - z4 * w) / s4;
    });

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_sigma") = mu_mu_mu_sigma,
        Named("mu_mu_sigma_sigma") = mu_mu_sigma_sigma,
        Named("mu_sigma_sigma_sigma") = mu_sigma_sigma_sigma,
        Named("sigma_sigma_sigma_sigma") = sigma_sigma_sigma_sigma
    );
}

// [[Rcpp::export]]
List gumbel_deriv4_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_sigma(n), mu_mu_sigma_sigma(n),
                  mu_sigma_sigma_sigma(n), sigma_sigma_sigma_sigma(n);
    bool sigma_is_scalar = (sigma.size() == 1);
    const double euler = 0.577215664901532860606512090082;
    double g[5];
    gumbel_gamma_moments(g);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double s4 = s * s * s * s;

        mu_mu_mu_mu[i] = -1.0 / s4;
        mu_mu_mu_sigma[i] = (4.0 - euler) / s4;
        mu_mu_sigma_sigma[i] = -(12.0 - 6.0 * euler + g[2]) / s4;
        mu_sigma_sigma_sigma[i] = (18.0 - 18.0 * euler + 9.0 * g[2] + g[3]) / s4;
        sigma_sigma_sigma_sigma[i] = (-18.0 - 36.0 * g[2] - 12.0 * g[3] - g[4]) / s4;
    });

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_sigma") = mu_mu_mu_sigma,
        Named("mu_mu_sigma_sigma") = mu_mu_sigma_sigma,
        Named("mu_sigma_sigma_sigma") = mu_sigma_sigma_sigma,
        Named("sigma_sigma_sigma_sigma") = sigma_sigma_sigma_sigma
    );
}
