#include <Rcpp.h>
using namespace Rcpp;

// Third/fourth-order derivatives of the Azzalini skew normal log-density.
//
// With z = (y - mu)/sigma and t = alpha z, the log-likelihood is
// l = const - log(sigma) - z^2/2 + log(Phi(t)). The derivatives of log(Phi)
// are polynomials in t and the inverse Mills ratio R = phi/Phi through
// R' = -R(t + R):
//   g1 = R
//   g2 = -R (t + R)
//   g3 = -g2 (t + R) - R (1 + g2)
//   g4 = -g3 (t + 2R) - 2 g2 (1 + g2)
// R is formed on the log scale, because below about t = -38 both phi and Phi
// underflow while the ratio stays finite and close to -t (same device as
// mills_ratio() on the R side).
//
// The (mu, sigma) ladder uses h_k, the k-th z-derivative of
// h(z) = -z^2/2 + log Phi(alpha z), with the location-scale generating rules
//   d/dmu    [A(z)/sigma^k] = -A'(z) / sigma^(k+1)
//   d/dsigma [A(z)/sigma^k] = -(z A'(z) + k A(z)) / sigma^(k+1);
// the alpha components differentiate g(t) directly through t = alpha z.

static inline void sn_g(double t, double g[5]) {
    double R = std::exp(R::dnorm4(t, 0.0, 1.0, 1) - R::pnorm5(t, 0.0, 1.0, 1, 1));
    double g2 = -R * (t + R);
    double g3 = -g2 * (t + R) - R * (1.0 + g2);
    double g4 = -g3 * (t + 2.0 * R) - 2.0 * g2 * (1.0 + g2);
    g[1] = R; g[2] = g2; g[3] = g3; g[4] = g4;
}

// [[Rcpp::export]]
List skewnormal_deriv3_cpp(NumericVector y, NumericVector mu,
                           NumericVector sigma, NumericVector alpha) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_sigma(n), mu_mu_alpha(n),
                  mu_sigma_sigma(n), mu_sigma_alpha(n), mu_alpha_alpha(n),
                  sigma_sigma_sigma(n), sigma_sigma_alpha(n),
                  sigma_alpha_alpha(n), alpha_alpha_alpha(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);
    bool alpha_is_scalar = (alpha.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double a = alpha_is_scalar ? alpha[0] : alpha[i];
        double s2 = s * s, s3 = s2 * s;
        double z = (y[i] - m) / s;
        double t = a * z;
        double z2 = z * z, z3 = z2 * z, t2 = t * t;

        double g[5];
        sn_g(t, g);
        double h1 = -z + a * g[1];
        double h2 = -1.0 + a * a * g[2];
        double h3 = a * a * a * g[3];

        mu_mu_mu[i] = -h3 / s3;
        mu_mu_sigma[i] = -(z * h3 + 2.0 * h2) / s3;
        mu_mu_alpha[i] = a * (2.0 * g[2] + t * g[3]) / s2;
        mu_sigma_sigma[i] = -(z2 * h3 + 4.0 * z * h2 + 2.0 * h1) / s3;
        mu_sigma_alpha[i] = (g[1] + 3.0 * t * g[2] + t2 * g[3]) / s2;
        mu_alpha_alpha[i] = -z * (2.0 * g[2] + t * g[3]) / s;
        sigma_sigma_sigma[i] = -(2.0 + 6.0 * z * h1 + 6.0 * z2 * h2 + z3 * h3) / s3;
        sigma_sigma_alpha[i] = z * (2.0 * g[1] + 4.0 * t * g[2] + t2 * g[3]) / s2;
        sigma_alpha_alpha[i] = -z2 * (2.0 * g[2] + t * g[3]) / s;
        alpha_alpha_alpha[i] = z3 * g[3];
    }

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_sigma") = mu_mu_sigma,
        Named("mu_mu_alpha") = mu_mu_alpha,
        Named("mu_sigma_sigma") = mu_sigma_sigma,
        Named("mu_sigma_alpha") = mu_sigma_alpha,
        Named("mu_alpha_alpha") = mu_alpha_alpha,
        Named("sigma_sigma_sigma") = sigma_sigma_sigma,
        Named("sigma_sigma_alpha") = sigma_sigma_alpha,
        Named("sigma_alpha_alpha") = sigma_alpha_alpha,
        Named("alpha_alpha_alpha") = alpha_alpha_alpha
    );
}

// [[Rcpp::export]]
List skewnormal_deriv4_cpp(NumericVector y, NumericVector mu,
                           NumericVector sigma, NumericVector alpha) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_sigma(n), mu_mu_mu_alpha(n),
                  mu_mu_sigma_sigma(n), mu_mu_sigma_alpha(n), mu_mu_alpha_alpha(n),
                  mu_sigma_sigma_sigma(n), mu_sigma_sigma_alpha(n),
                  mu_sigma_alpha_alpha(n), mu_alpha_alpha_alpha(n),
                  sigma_sigma_sigma_sigma(n), sigma_sigma_sigma_alpha(n),
                  sigma_sigma_alpha_alpha(n), sigma_alpha_alpha_alpha(n),
                  alpha_alpha_alpha_alpha(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);
    bool alpha_is_scalar = (alpha.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double a = alpha_is_scalar ? alpha[0] : alpha[i];
        double s2 = s * s, s3 = s2 * s, s4 = s2 * s2;
        double z = (y[i] - m) / s;
        double t = a * z;
        double z2 = z * z, z3 = z2 * z, z4 = z2 * z2;
        double t2 = t * t, t3 = t2 * t;

        double g[5];
        sn_g(t, g);
        double h1 = -z + a * g[1];
        double h2 = -1.0 + a * a * g[2];
        double h3 = a * a * a * g[3];
        double h4 = a * a * a * a * g[4];

        mu_mu_mu_mu[i] = h4 / s4;
        mu_mu_mu_sigma[i] = (z * h4 + 3.0 * h3) / s4;
        mu_mu_mu_alpha[i] = -a * a * (3.0 * g[3] + t * g[4]) / s3;
        mu_mu_sigma_sigma[i] = (z2 * h4 + 6.0 * z * h3 + 6.0 * h2) / s4;
        mu_mu_sigma_alpha[i] = -a * (4.0 * g[2] + 5.0 * t * g[3] + t2 * g[4]) / s3;
        mu_mu_alpha_alpha[i] = (2.0 * g[2] + 4.0 * t * g[3] + t2 * g[4]) / s2;
        mu_sigma_sigma_sigma[i] = (z3 * h4 + 9.0 * z2 * h3 + 18.0 * z * h2
                                   + 6.0 * h1) / s4;
        mu_sigma_sigma_alpha[i] = -(2.0 * g[1] + 10.0 * t * g[2]
                                    + 7.0 * t2 * g[3] + t3 * g[4]) / s3;
        mu_sigma_alpha_alpha[i] = z * (4.0 * g[2] + 5.0 * t * g[3] + t2 * g[4]) / s2;
        mu_alpha_alpha_alpha[i] = -z2 * (3.0 * g[3] + t * g[4]) / s;
        sigma_sigma_sigma_sigma[i] = (6.0 + 24.0 * z * h1 + 36.0 * z2 * h2
                                      + 12.0 * z3 * h3 + z4 * h4) / s4;
        sigma_sigma_sigma_alpha[i] = -z * (6.0 * g[1] + 18.0 * t * g[2]
                                           + 9.0 * t2 * g[3] + t3 * g[4]) / s3;
        sigma_sigma_alpha_alpha[i] = z2 * (6.0 * g[2] + 6.0 * t * g[3]
                                           + t2 * g[4]) / s2;
        sigma_alpha_alpha_alpha[i] = -z3 * (3.0 * g[3] + t * g[4]) / s;
        alpha_alpha_alpha_alpha[i] = z4 * g[4];
    }

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_sigma") = mu_mu_mu_sigma,
        Named("mu_mu_mu_alpha") = mu_mu_mu_alpha,
        Named("mu_mu_sigma_sigma") = mu_mu_sigma_sigma,
        Named("mu_mu_sigma_alpha") = mu_mu_sigma_alpha,
        Named("mu_mu_alpha_alpha") = mu_mu_alpha_alpha,
        Named("mu_sigma_sigma_sigma") = mu_sigma_sigma_sigma,
        Named("mu_sigma_sigma_alpha") = mu_sigma_sigma_alpha,
        Named("mu_sigma_alpha_alpha") = mu_sigma_alpha_alpha,
        Named("mu_alpha_alpha_alpha") = mu_alpha_alpha_alpha,
        Named("sigma_sigma_sigma_sigma") = sigma_sigma_sigma_sigma,
        Named("sigma_sigma_sigma_alpha") = sigma_sigma_sigma_alpha,
        Named("sigma_sigma_alpha_alpha") = sigma_sigma_alpha_alpha,
        Named("sigma_alpha_alpha_alpha") = sigma_alpha_alpha_alpha,
        Named("alpha_alpha_alpha_alpha") = alpha_alpha_alpha_alpha
    );
}
