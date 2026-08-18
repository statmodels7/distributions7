#include <Rcpp.h>
#include "d7_par.h"
using namespace Rcpp;

// Third/fourth-order derivatives of the Negative Binomial (NB2) log-mass,
// derived from the gradient/Hessian by repeated differentiation and validated
// numerically. s = theta + mu.
//
// E[psigamma(Y + theta, deriv)] for the expected pure-theta derivatives, summed
// over the support up to a far-tail quantile with the same pmf recurrence used
// by the expected Hessian; the residual tail mass (~1e-12) is negligible.
static double nb_E_psigamma(double mu, double theta, double deriv) {
    double ratio = mu / (theta + mu);
    double kq = R::qnbinom_mu(1.0 - 1e-12, theta, mu, 1, 0);
    int kmax = (int) std::max(100.0, kq) + 1;

    double s = 0.0, cum = 0.0;
    double log_p0 = theta * (std::log(theta) - std::log(theta + mu));

    if (log_p0 > -700.0) {
        double pk = std::exp(log_p0);
        for (int k = 0; k <= kmax; ++k) {
            s += R::psigamma(k + theta, deriv) * pk;
            cum += pk;
            pk *= (k + theta) / (k + 1.0) * ratio;
        }
    } else {
        for (int k = 0; k <= kmax; ++k) {
            double pk = R::dnbinom_mu(k, theta, mu, 0);
            s += R::psigamma(k + theta, deriv) * pk;
            cum += pk;
        }
    }
    if (cum < 1.0) s += R::psigamma(kmax + 1 + theta, deriv) * (1.0 - cum);
    return s;
}

// [[Rcpp::export]]
List negbin_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector theta,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_theta(n), mu_theta_theta(n), theta_theta_theta(n);
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_s ? mu[0] : mu[i];
        double th = th_s ? theta[0] : theta[i];
        double yi = y[i];
        double s = th + m, s2 = s * s, s3 = s2 * s;
        double m2 = m * m, m3 = m2 * m, th2 = th * th;

        mu_mu_mu[i] = -2.0 * (yi + th) / s3 + 2.0 * yi / m3;
        mu_mu_theta[i] = 1.0 / s2 - 2.0 * (yi + th) / s3;
        mu_theta_theta[i] = -2.0 * (yi - m) / s3;
        theta_theta_theta[i] = R::psigamma(yi + th, 2.0) - R::psigamma(th, 2.0)
            - m * (2.0 * th + m) / (th2 * s2) - 2.0 * (yi - m) / s3;
    });

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_theta") = mu_mu_theta,
        Named("mu_theta_theta") = mu_theta_theta,
        Named("theta_theta_theta") = theta_theta_theta
    );
}

// [[Rcpp::export]]
List negbin_deriv3_expected_cpp(NumericVector y, NumericVector mu, NumericVector theta,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_theta(n), mu_theta_theta(n), theta_theta_theta(n);
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);
    bool both_scalar = mu_s && th_s;

    // the scalar-case constant lives OUT here; the per-iteration copy is
    // LOCAL to the lambda, or two threads would race on it
    double ttt0 = 0;
    if (both_scalar) {
        double m = mu[0], th = theta[0], s = th + m;
        ttt0 = nb_E_psigamma(m, th, 2.0) - R::psigamma(th, 2.0) - m * (2.0 * th + m) / (th * th * s * s);
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_s ? mu[0] : mu[i];
        double th = th_s ? theta[0] : theta[i];
        double s = th + m, s2 = s * s;
        double ttt = ttt0;
        if (!both_scalar) {
            ttt = nb_E_psigamma(m, th, 2.0) - R::psigamma(th, 2.0) - m * (2.0 * th + m) / (th * th * s2);
        }
        mu_mu_mu[i] = -2.0 / s2 + 2.0 / (m * m);
        mu_mu_theta[i] = -1.0 / s2;
        mu_theta_theta[i] = 0.0;
        theta_theta_theta[i] = ttt;
    });

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_theta") = mu_mu_theta,
        Named("mu_theta_theta") = mu_theta_theta,
        Named("theta_theta_theta") = theta_theta_theta
    );
}

// [[Rcpp::export]]
List negbin_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector theta,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_theta(n), mu_mu_theta_theta(n),
                  mu_theta_theta_theta(n), theta_theta_theta_theta(n);
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_s ? mu[0] : mu[i];
        double th = th_s ? theta[0] : theta[i];
        double yi = y[i];
        double s = th + m, s2 = s * s, s3 = s2 * s, s4 = s2 * s2;
        double m2 = m * m, m4 = m2 * m2, th3 = th * th * th;

        mu_mu_mu_mu[i] = 6.0 * (yi + th) / s4 - 6.0 * yi / m4;
        mu_mu_mu_theta[i] = -2.0 / s3 + 6.0 * (yi + th) / s4;
        mu_mu_theta_theta[i] = -4.0 / s3 + 6.0 * (yi + th) / s4;
        mu_theta_theta_theta[i] = 6.0 * (yi - m) / s4;
        theta_theta_theta_theta[i] = R::psigamma(yi + th, 3.0) - R::psigamma(th, 3.0)
            - 2.0 * m * (th * s - (2.0 * th + m) * (2.0 * th + m)) / (th3 * s3)
            + 6.0 * (yi - m) / s4;
    });

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_theta") = mu_mu_mu_theta,
        Named("mu_mu_theta_theta") = mu_mu_theta_theta,
        Named("mu_theta_theta_theta") = mu_theta_theta_theta,
        Named("theta_theta_theta_theta") = theta_theta_theta_theta
    );
}

// [[Rcpp::export]]
List negbin_deriv4_expected_cpp(NumericVector y, NumericVector mu, NumericVector theta,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_theta(n), mu_mu_theta_theta(n),
                  mu_theta_theta_theta(n), theta_theta_theta_theta(n);
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);
    bool both_scalar = mu_s && th_s;

    double tttt0 = 0;
    if (both_scalar) {
        double m = mu[0], th = theta[0], s = th + m;
        double th3 = th * th * th, s3 = s * s * s;
        tttt0 = nb_E_psigamma(m, th, 3.0) - R::psigamma(th, 3.0)
            - 2.0 * m * (th * s - (2.0 * th + m) * (2.0 * th + m)) / (th3 * s3);
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_s ? mu[0] : mu[i];
        double th = th_s ? theta[0] : theta[i];
        double s = th + m, s3 = s * s * s;
        double tttt = tttt0;
        if (!both_scalar) {
            double th3 = th * th * th;
            tttt = nb_E_psigamma(m, th, 3.0) - R::psigamma(th, 3.0)
                - 2.0 * m * (th * s - (2.0 * th + m) * (2.0 * th + m)) / (th3 * s3);
        }
        mu_mu_mu_mu[i] = 6.0 / s3 - 6.0 / (m * m * m);
        mu_mu_mu_theta[i] = 4.0 / s3;
        mu_mu_theta_theta[i] = 2.0 / s3;
        mu_theta_theta_theta[i] = 0.0;
        theta_theta_theta_theta[i] = tttt;
    });

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_theta") = mu_mu_mu_theta,
        Named("mu_mu_theta_theta") = mu_mu_theta_theta,
        Named("mu_theta_theta_theta") = mu_theta_theta_theta,
        Named("theta_theta_theta_theta") = theta_theta_theta_theta
    );
}
