#include <Rcpp.h>
#include "d7_par.h"
using namespace Rcpp;

// Third/fourth-order derivatives of the Inverse-Gaussian log-density (mean-
// dispersion parameterization), transcribed from the Wolfram output.

// [[Rcpp::export]]
List invgauss_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector phi,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_phi(n), mu_phi_phi(n), phi_phi_phi(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool phi_is_scalar = (phi.size() == 1);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double p = phi_is_scalar ? phi[0] : phi[i];
        double m2 = m * m, m3 = m2 * m, m4 = m2 * m2, m5 = m4 * m;
        double p2 = p * p, p3 = p2 * p, p4 = p2 * p2;
        double r = m - y[i];

        mu_mu_mu[i] = -6.0 * (m - 2.0 * y[i]) / (m5 * p);
        mu_mu_phi[i] = (-2.0 * m + 3.0 * y[i]) / (m4 * p2);
        mu_phi_phi[i] = -2.0 * r / (m3 * p3);
        phi_phi_phi[i] = (-p + 3.0 * r * r / (m2 * y[i])) / p4;
    });

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_phi") = mu_mu_phi,
        Named("mu_phi_phi") = mu_phi_phi,
        Named("phi_phi_phi") = phi_phi_phi
    );
}

// [[Rcpp::export]]
List invgauss_deriv3_expected_cpp(NumericVector y, NumericVector mu, NumericVector phi,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_phi(n), mu_phi_phi(n), phi_phi_phi(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool phi_is_scalar = (phi.size() == 1);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double p = phi_is_scalar ? phi[0] : phi[i];
        double m3 = m * m * m, m4 = m3 * m;
        double p2 = p * p, p3 = p2 * p;

        mu_mu_mu[i] = 6.0 / (m4 * p);
        mu_mu_phi[i] = 1.0 / (m3 * p2);
        mu_phi_phi[i] = 0.0;
        phi_phi_phi[i] = 2.0 / p3;
    });

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_phi") = mu_mu_phi,
        Named("mu_phi_phi") = mu_phi_phi,
        Named("phi_phi_phi") = phi_phi_phi
    );
}

// [[Rcpp::export]]
List invgauss_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector phi,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_phi(n), mu_mu_phi_phi(n),
                  mu_phi_phi_phi(n), phi_phi_phi_phi(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool phi_is_scalar = (phi.size() == 1);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double p = phi_is_scalar ? phi[0] : phi[i];
        double m2 = m * m, m3 = m2 * m, m4 = m2 * m2, m5 = m4 * m, m6 = m4 * m2;
        double p2 = p * p, p3 = p2 * p, p4 = p2 * p2, p5 = p4 * p;
        double r = m - y[i];

        mu_mu_mu_mu[i] = (24.0 * m - 60.0 * y[i]) / (m6 * p);
        mu_mu_mu_phi[i] = 6.0 * (m - 2.0 * y[i]) / (m5 * p2);
        mu_mu_phi_phi[i] = (4.0 * m - 6.0 * y[i]) / (m4 * p3);
        mu_phi_phi_phi[i] = 6.0 * r / (m3 * p4);
        phi_phi_phi_phi[i] = 3.0 * (p - 4.0 * r * r / (m2 * y[i])) / p5;
    });

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_phi") = mu_mu_mu_phi,
        Named("mu_mu_phi_phi") = mu_mu_phi_phi,
        Named("mu_phi_phi_phi") = mu_phi_phi_phi,
        Named("phi_phi_phi_phi") = phi_phi_phi_phi
    );
}

// [[Rcpp::export]]
List invgauss_deriv4_expected_cpp(NumericVector y, NumericVector mu, NumericVector phi,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_phi(n), mu_mu_phi_phi(n),
                  mu_phi_phi_phi(n), phi_phi_phi_phi(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool phi_is_scalar = (phi.size() == 1);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double p = phi_is_scalar ? phi[0] : phi[i];
        double m3 = m * m * m, m4 = m3 * m, m5 = m4 * m;
        double p2 = p * p, p3 = p2 * p, p4 = p2 * p2;

        mu_mu_mu_mu[i] = -36.0 / (m5 * p);
        mu_mu_mu_phi[i] = -6.0 / (m4 * p2);
        mu_mu_phi_phi[i] = -2.0 / (m3 * p3);
        mu_phi_phi_phi[i] = 0.0;
        phi_phi_phi_phi[i] = -9.0 / p4;
    });

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_phi") = mu_mu_mu_phi,
        Named("mu_mu_phi_phi") = mu_mu_phi_phi,
        Named("mu_phi_phi_phi") = mu_phi_phi_phi,
        Named("phi_phi_phi_phi") = phi_phi_phi_phi
    );
}
