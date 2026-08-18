#include <Rcpp.h>
#include "d7_par.h"
using namespace Rcpp;

// E[trigamma(Y + theta)] for Y ~ NB2(mu, theta), needed by the expected hessian.
// Sums trigamma(k + theta) * P(Y = k) up to a far-tail quantile, using the
// pmf recurrence p_{k+1} = p_k * (k + theta) / (k + 1) * mu / (theta + mu),
// with a log-space fallback when p_0 underflows.
static double nb_E_trigamma(double mu, double theta) {
    double ratio = mu / (theta + mu);
    double kq = R::qnbinom_mu(1.0 - 1e-12, theta, mu, 1, 0);
    int kmax = (int) std::max(100.0, kq) + 1;

    double s = 0.0, cum = 0.0;
    double log_p0 = theta * (std::log(theta) - std::log(theta + mu));

    if (log_p0 > -700.0) {
        double pk = std::exp(log_p0);
        for (int k = 0; k <= kmax; ++k) {
            s += R::trigamma(k + theta) * pk;
            cum += pk;
            pk *= (k + theta) / (k + 1.0) * ratio;
        }
    } else {
        for (int k = 0; k <= kmax; ++k) {
            double pk = R::dnbinom_mu(k, theta, mu, 0);
            s += R::trigamma(k + theta) * pk;
            cum += pk;
        }
    }

    // Tail bound: trigamma is decreasing, so the missing mass contributes at
    // most trigamma(kmax + 1 + theta) * (1 - cum).
    if (cum < 1.0) s += R::trigamma(kmax + 1 + theta) * (1.0 - cum);
    return s;
}

// [[Rcpp::export]]
List negbin_gradient_cpp(NumericVector y, NumericVector mu, NumericVector theta,
                        int threads = 1) {
    int n = y.size();
    NumericVector grad_mu(n);
    NumericVector grad_theta(n);

    bool mu_is_scalar = (mu.size() == 1);
    bool theta_is_scalar = (theta.size() == 1);
    bool both_scalar = mu_is_scalar && theta_is_scalar;

    // the scalar-case constants live OUT here; the per-iteration copies are
    // LOCAL to the lambda, or two threads would race on them
    double m0 = 0, th0 = 0, th_plus_mu0 = 0, digamma_th0 = 0, log_frac0 = 0;

    if (both_scalar) {
        m0 = mu[0];
        th0 = theta[0];
        th_plus_mu0 = th0 + m0;
        digamma_th0 = R::digamma(th0);
        log_frac0 = std::log(th0 / th_plus_mu0);
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m0, th = th0, th_plus_mu = th_plus_mu0,
               digamma_th = digamma_th0, log_frac = log_frac0;
        if (!both_scalar) {
            m = mu_is_scalar ? mu[0] : mu[i];
            th = theta_is_scalar ? theta[0] : theta[i];
            th_plus_mu = th + m;
            digamma_th = R::digamma(th);
            log_frac = std::log(th / th_plus_mu);
        }

        grad_mu[i] = (th / th_plus_mu) * (y[i] / m - 1.0);
        grad_theta[i] = R::digamma(y[i] + th) - digamma_th + log_frac + (m - y[i]) / th_plus_mu;
    });

    return List::create(Named("mu") = grad_mu, Named("theta") = grad_theta);
}

// [[Rcpp::export]]
List negbin_hessian_cpp(NumericVector y, NumericVector mu, NumericVector theta,
                        int threads = 1) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    NumericVector hess_theta_theta(n);
    NumericVector hess_mu_theta(n);

    bool mu_is_scalar = (mu.size() == 1);
    bool theta_is_scalar = (theta.size() == 1);
    bool both_scalar = mu_is_scalar && theta_is_scalar;

    double m0 = 0, th0 = 0, th_plus_mu0 = 0, th_plus_mu20 = 0,
           trigamma_th0 = 0, mid0 = 0;

    if (both_scalar) {
        m0 = mu[0];
        th0 = theta[0];
        th_plus_mu0 = th0 + m0;
        th_plus_mu20 = th_plus_mu0 * th_plus_mu0;
        trigamma_th0 = R::trigamma(th0);
        mid0 = m0 / (th0 * th_plus_mu0);
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m0, th = th0, th_plus_mu = th_plus_mu0,
               th_plus_mu2 = th_plus_mu20, trigamma_th = trigamma_th0,
               mid = mid0;
        if (!both_scalar) {
            m = mu_is_scalar ? mu[0] : mu[i];
            th = theta_is_scalar ? theta[0] : theta[i];
            th_plus_mu = th + m;
            th_plus_mu2 = th_plus_mu * th_plus_mu;
            trigamma_th = R::trigamma(th);
            mid = m / (th * th_plus_mu);
        }

        double res = y[i] - m;

        hess_mu_mu[i] = (y[i] + th) / th_plus_mu2 - y[i] / (m * m);
        hess_theta_theta[i] = R::trigamma(y[i] + th) - trigamma_th + mid + res / th_plus_mu2;
        hess_mu_theta[i] = res / th_plus_mu2;
    });

    return List::create(
        Named("mu_mu") = hess_mu_mu,
        Named("theta_theta") = hess_theta_theta,
        Named("mu_theta") = hess_mu_theta
    );
}

// [[Rcpp::export]]
List negbin_expected_hessian_cpp(NumericVector y, NumericVector mu, NumericVector theta,
                        int threads = 1) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    NumericVector hess_theta_theta(n);
    NumericVector hess_mu_theta(n);

    bool mu_is_scalar = (mu.size() == 1);
    bool theta_is_scalar = (theta.size() == 1);
    bool both_scalar = mu_is_scalar && theta_is_scalar;

    double hmm0 = 0, htt0 = 0;

    if (both_scalar) {
        double m = mu[0];
        double th = theta[0];
        double th_plus_mu = th + m;
        hmm0 = -th / (m * th_plus_mu);
        htt0 = nb_E_trigamma(m, th) - R::trigamma(th) + m / (th * th_plus_mu);
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double hmm = hmm0, htt = htt0;
        if (!both_scalar) {
            double m = mu_is_scalar ? mu[0] : mu[i];
            double th = theta_is_scalar ? theta[0] : theta[i];
            double th_plus_mu = th + m;
            hmm = -th / (m * th_plus_mu);
            htt = nb_E_trigamma(m, th) - R::trigamma(th) + m / (th * th_plus_mu);
        }

        hess_mu_mu[i] = hmm;
        hess_theta_theta[i] = htt;
        hess_mu_theta[i] = 0.0;
    });

    return List::create(
        Named("mu_mu") = hess_mu_mu,
        Named("theta_theta") = hess_theta_theta,
        Named("mu_theta") = hess_mu_theta
    );
}
