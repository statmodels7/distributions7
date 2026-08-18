#include <Rcpp.h>
#include "d7_par.h"
using namespace Rcpp;

// Third/fourth-order derivatives of the Beta log-density (mean/precision
// parameterization, a = mu*phi, b = (1-mu)*phi), transcribed from the Wolfram
// output. These derivatives do not depend on y (the y-terms of the log-density
// are linear in the parameters), so the observed and expected derivatives
// coincide -- a single kernel serves both.

// [[Rcpp::export]]
List beta_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector phi,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_phi(n), mu_phi_phi(n), phi_phi_phi(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool phi_is_scalar = (phi.size() == 1);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double p = phi_is_scalar ? phi[0] : phi[i];
        double a = m * p, b = p - m * p;
        double p2 = p * p, p3 = p2 * p;
        double Ta = R::trigamma(a), Tb = R::trigamma(b);
        double P2a = R::psigamma(a, 2.0), P2b = R::psigamma(b, 2.0);
        double P2phi = R::psigamma(p, 2.0);

        mu_mu_mu[i] = p3 * (-P2a + P2b);
        mu_mu_phi[i] = p * (-(m * p * P2a) + (m - 1.0) * p * P2b - 2.0 * (Ta + Tb));
        mu_phi_phi[i] = -(m * m * p * P2a) + (m - 1.0) * (m - 1.0) * p * P2b - 2.0 * m * Ta - 2.0 * (m - 1.0) * Tb;
        phi_phi_phi[i] = P2phi - m * m * m * P2a + (m - 1.0) * (m - 1.0) * (m - 1.0) * P2b;
    });

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_phi") = mu_mu_phi,
        Named("mu_phi_phi") = mu_phi_phi,
        Named("phi_phi_phi") = phi_phi_phi
    );
}

// [[Rcpp::export]]
List beta_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector phi,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_phi(n), mu_mu_phi_phi(n),
                  mu_phi_phi_phi(n), phi_phi_phi_phi(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool phi_is_scalar = (phi.size() == 1);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double p = phi_is_scalar ? phi[0] : phi[i];
        double a = m * p, b = p - m * p;
        double om = m - 1.0;
        double p2 = p * p, p4 = p2 * p2;
        double Ta = R::trigamma(a), Tb = R::trigamma(b);
        double P2a = R::psigamma(a, 2.0), P2b = R::psigamma(b, 2.0);
        double P3a = R::psigamma(a, 3.0), P3b = R::psigamma(b, 3.0), P3phi = R::psigamma(p, 3.0);

        mu_mu_mu_mu[i] = -(p4 * (P3a + P3b));
        mu_mu_mu_phi[i] = p2 * (-3.0 * P2a - m * p * P3a + 3.0 * P2b - om * p * P3b);
        mu_mu_phi_phi[i] = p * (-4.0 * m * P2a - m * m * p * P3a + om * (4.0 * P2b - om * p * P3b)) - 2.0 * (Ta + Tb);
        mu_phi_phi_phi[i] = -3.0 * m * m * P2a - m * m * m * p * P3a + om * om * (3.0 * P2b - om * p * P3b);
        phi_phi_phi_phi[i] = P3phi - m * m * m * m * P3a - om * om * om * om * P3b;
    });

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_phi") = mu_mu_mu_phi,
        Named("mu_mu_phi_phi") = mu_mu_phi_phi,
        Named("mu_phi_phi_phi") = mu_phi_phi_phi,
        Named("phi_phi_phi_phi") = phi_phi_phi_phi
    );
}
