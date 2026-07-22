#include <Rcpp.h>
using namespace Rcpp;

// Third/fourth-order observed derivatives of the Logistic log-density.
//
// The Wolfram output for this model is written directly in exp((mu - y)/sigma),
// with numerators and denominators up to the fourth power of it: transcribed
// literally it overflows for moderately large standardized residuals. The
// expressions below are an equivalent rewriting in the sigmoid, which keeps
// every intermediate quantity inside [0, 1].
//
// With z = (y - mu)/sigma, t = 1/(1 + exp(-z)) and u = 1 - t, the log-density is
// l = -log(sigma) + g(z) with g(z) = -z - 2*log(1 + exp(-z)), whose derivatives
// with respect to z are
//
//   g1 = 1 - 2t,  g2 = -2tu,  g3 = -2tu(1 - 2t),  g4 = -2tu(1 - 6tu).
//
// Because z = (y - mu)/sigma, differentiating a term of the form A(z)/sigma^k
// obeys the two rules
//
//   d/dmu    [A/sigma^k] = -A'/sigma^(k+1)
//   d/dsigma [A/sigma^k] = -(z*A' + k*A)/sigma^(k+1)
//
// and applying them repeatedly gives every component as a polynomial in z whose
// coefficients are the g's. The mixed partials agree whichever order the two
// rules are applied in, which is what pins the coefficients below.

// Sigmoid evaluated on the stable branch: the exponential never sees a positive
// argument, so no overflow is possible whatever the residual.
static inline void sigmoid_pair(double z, double& t, double& u) {
    if (z >= 0.0) {
        double e = std::exp(-z);
        t = 1.0 / (1.0 + e);
        u = e / (1.0 + e);
    } else {
        double e = std::exp(z);
        t = e / (1.0 + e);
        u = 1.0 / (1.0 + e);
    }
}

// [[Rcpp::export]]
List logistic_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector sigma) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_sigma(n), mu_sigma_sigma(n), sigma_sigma_sigma(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double s3 = s * s * s;

        double z = (y[i] - m) / s;
        double z2 = z * z, z3 = z2 * z;

        double t, u;
        sigmoid_pair(z, t, u);
        double tu = t * u;

        double g1 = u - t;
        double g2 = -2.0 * tu;
        double g3 = g2 * g1;

        mu_mu_mu[i] = -g3 / s3;
        mu_mu_sigma[i] = -(2.0 * g2 + z * g3) / s3;
        mu_sigma_sigma[i] = -(2.0 * g1 + 4.0 * z * g2 + z2 * g3) / s3;
        sigma_sigma_sigma[i] = -(2.0 + 6.0 * z * g1 + 6.0 * z2 * g2 + z3 * g3) / s3;
    }

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_sigma") = mu_mu_sigma,
        Named("mu_sigma_sigma") = mu_sigma_sigma,
        Named("sigma_sigma_sigma") = sigma_sigma_sigma
    );
}

// [[Rcpp::export]]
List logistic_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector sigma) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_sigma(n), mu_mu_sigma_sigma(n),
                  mu_sigma_sigma_sigma(n), sigma_sigma_sigma_sigma(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double s2 = s * s, s4 = s2 * s2;

        double z = (y[i] - m) / s;
        double z2 = z * z, z3 = z2 * z, z4 = z2 * z2;

        double t, u;
        sigmoid_pair(z, t, u);
        double tu = t * u;

        double g1 = u - t;
        double g2 = -2.0 * tu;
        double g3 = g2 * g1;
        double g4 = g2 * (1.0 - 6.0 * tu);

        mu_mu_mu_mu[i] = g4 / s4;
        mu_mu_mu_sigma[i] = (3.0 * g3 + z * g4) / s4;
        mu_mu_sigma_sigma[i] = (6.0 * g2 + 6.0 * z * g3 + z2 * g4) / s4;
        mu_sigma_sigma_sigma[i] = (6.0 * g1 + 18.0 * z * g2 + 9.0 * z2 * g3 + z3 * g4) / s4;
        sigma_sigma_sigma_sigma[i] =
            (6.0 + 24.0 * z * g1 + 36.0 * z2 * g2 + 12.0 * z3 * g3 + z4 * g4) / s4;
    }

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_sigma") = mu_mu_mu_sigma,
        Named("mu_mu_sigma_sigma") = mu_mu_sigma_sigma,
        Named("mu_sigma_sigma_sigma") = mu_sigma_sigma_sigma,
        Named("sigma_sigma_sigma_sigma") = sigma_sigma_sigma_sigma
    );
}
