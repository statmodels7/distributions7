#include <Rcpp.h>
using namespace Rcpp;

// Third/fourth-order derivatives of the Weibull log-density in the gamlss WEI
// parametrization (mu = scale, sigma = shape). With u = (y/mu)^sigma and
// L = log(y/mu), du/dmu = -sigma*u/mu and du/dsigma = u*L, so every derivative
// stays a polynomial in u and L*u over powers of mu and sigma; the ladders are
// derived by hand from the gradient the package already ships.

// Derivatives of Gamma at 2, i.e. E[u (log u)^k] for a standard exponential u,
// assembled from polygamma values at 2 by the moment-cumulant relations.
static void weibull_gamma_moments(double g[5]) {
    const double euler = 0.577215664901532860606512090082;
    const double zeta3 = 1.202056903159594285399738161511;
    const double psi  = 1.0 - euler;                       // psi(2)
    const double psi1 = M_PI * M_PI / 6.0 - 1.0;           // psi'(2)
    const double psi2 = 2.0 - 2.0 * zeta3;                 // psi''(2)
    const double psi3 = M_PI * M_PI * M_PI * M_PI / 15.0 - 6.0;  // psi'''(2)
    g[0] = 1.0;
    g[1] = psi;
    g[2] = psi * psi + psi1;
    g[3] = psi * psi * psi + 3.0 * psi * psi1 + psi2;
    g[4] = psi * psi * psi * psi + 6.0 * psi * psi * psi1
         + 4.0 * psi * psi2 + 3.0 * psi1 * psi1 + psi3;
}

// [[Rcpp::export]]
List weibull_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector sigma) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_sigma(n), mu_sigma_sigma(n), sigma_sigma_sigma(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double m2 = m * m, m3 = m2 * m, s3 = s * s * s;
        double L = std::log(y[i] / m);
        double u = std::exp(s * L);

        mu_mu_mu[i] = s * (-2.0 + (1.0 + s) * (2.0 + s) * u) / m3;
        mu_mu_sigma[i] = (1.0 - (1.0 + 2.0 * s) * u - s * (1.0 + s) * L * u) / m2;
        mu_sigma_sigma[i] = L * u * (2.0 + s * L) / m;
        sigma_sigma_sigma[i] = 2.0 / s3 - u * L * L * L;
    }

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_sigma") = mu_mu_sigma,
        Named("mu_sigma_sigma") = mu_sigma_sigma,
        Named("sigma_sigma_sigma") = sigma_sigma_sigma
    );
}

// [[Rcpp::export]]
List weibull_deriv3_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_sigma(n), mu_sigma_sigma(n), sigma_sigma_sigma(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);
    double g[5];
    weibull_gamma_moments(g);

    // E[u] = 1 and E[u L^k] = Gamma^(k)(2) / sigma^k, since u is standard
    // exponential whatever the parameters and L = log(u)/sigma.
    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double m2 = m * m, m3 = m2 * m, s3 = s * s * s;

        mu_mu_mu[i] = s * s * (s + 3.0) / m3;
        mu_mu_sigma[i] = (-2.0 * s - (1.0 + s) * g[1]) / m2;
        mu_sigma_sigma[i] = (2.0 * g[1] + g[2]) / (s * m);
        sigma_sigma_sigma[i] = (2.0 - g[3]) / s3;
    }

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_sigma") = mu_mu_sigma,
        Named("mu_sigma_sigma") = mu_sigma_sigma,
        Named("sigma_sigma_sigma") = sigma_sigma_sigma
    );
}

// [[Rcpp::export]]
List weibull_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector sigma) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_sigma(n), mu_mu_sigma_sigma(n),
                  mu_sigma_sigma_sigma(n), sigma_sigma_sigma_sigma(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double m2 = m * m, m3 = m2 * m, m4 = m2 * m2, s4 = s * s * s * s;
        double L = std::log(y[i] / m);
        double u = std::exp(s * L);
        double L2 = L * L;

        mu_mu_mu_mu[i] = s * (6.0 - (1.0 + s) * (2.0 + s) * (3.0 + s) * u) / m4;
        mu_mu_mu_sigma[i] = (-2.0 + (2.0 + 6.0 * s + 3.0 * s * s) * u
                             + s * (1.0 + s) * (2.0 + s) * L * u) / m3;
        mu_mu_sigma_sigma[i] = -u * (2.0 + 2.0 * (1.0 + 2.0 * s) * L
                                     + s * (1.0 + s) * L2) / m2;
        mu_sigma_sigma_sigma[i] = L2 * u * (3.0 + s * L) / m;
        sigma_sigma_sigma_sigma[i] = -6.0 / s4 - u * L2 * L2;
    }

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_sigma") = mu_mu_mu_sigma,
        Named("mu_mu_sigma_sigma") = mu_mu_sigma_sigma,
        Named("mu_sigma_sigma_sigma") = mu_sigma_sigma_sigma,
        Named("sigma_sigma_sigma_sigma") = sigma_sigma_sigma_sigma
    );
}

// [[Rcpp::export]]
List weibull_deriv4_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_sigma(n), mu_mu_sigma_sigma(n),
                  mu_sigma_sigma_sigma(n), sigma_sigma_sigma_sigma(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);
    double g[5];
    weibull_gamma_moments(g);

    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double m2 = m * m, m3 = m2 * m, m4 = m2 * m2, s2 = s * s, s4 = s2 * s2;

        mu_mu_mu_mu[i] = s * (6.0 - (1.0 + s) * (2.0 + s) * (3.0 + s)) / m4;
        mu_mu_mu_sigma[i] = (6.0 * s + 3.0 * s2
                             + (1.0 + s) * (2.0 + s) * g[1]) / m3;
        mu_mu_sigma_sigma[i] = -(2.0 + 2.0 * (1.0 + 2.0 * s) * g[1] / s
                                 + (1.0 + s) * g[2] / s) / m2;
        mu_sigma_sigma_sigma[i] = (3.0 * g[2] + g[3]) / (s2 * m);
        sigma_sigma_sigma_sigma[i] = (-6.0 - g[4]) / s4;
    }

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_sigma") = mu_mu_mu_sigma,
        Named("mu_mu_sigma_sigma") = mu_mu_sigma_sigma,
        Named("mu_sigma_sigma_sigma") = mu_sigma_sigma_sigma,
        Named("sigma_sigma_sigma_sigma") = sigma_sigma_sigma_sigma
    );
}
