#include <Rcpp.h>
#include "d7_par.h"
#include "psi_diff.h"
using namespace Rcpp;

// THE DISPERSION AT LARGE theta.
//
// As theta grows the negative binomial tends to the Poisson and every
// derivative in theta vanishes, so each is written as a sum of terms that
// cancel to leading order.  Measured, the score's four terms cancel PAIRWISE:
// psi(y+th) - psi(th) is y/th, log(th/(th+mu)) is -mu/th and (mu-y)/(th+mu) is
// (mu-y)/th, and the three sum to zero, so the value is O(1/th^2) computed
// from terms of size 1/th -- a cancellation of order theta, not of order
// theta/y as a look at the digamma difference alone suggests.  The direct form
// is wrong by 1.0e-03 at theta = 1e6, by 4.4 at 1e7 and CHANGES SIGN at 1e8.
//
// And a fit reaches there routinely: on 2000 counts with mu = 4 drawn at a
// true theta of 100, `fit_distrib` reports 1.6e+07; on Poisson counts it
// reports 2.3e+05.  Where it stops in that limit is therefore decided by
// which wrong value happens to cross the tolerance.
//
// Both quantities are rewritten so that each cancellation is performed
// SYMBOLICALLY and what is left is evaluated directly.  With a = theta,
// b = theta + y and c = theta + mu:
//
//   dl/dtheta   = [psi(b) - psi(a) - log1p(y/a)] + [log1p(w) - w],
//                 w = (y - mu)/c
//   d2l/dtheta2 = (y - mu)^2/(b c^2)
//                 + [psi'(b) - psi'(a) + y/(a b)]
//
// The first bracket of the score is y/(2ab) + ... and the second is -w^2/2 +
// ..., each with its own series; the leading three terms of the Hessian
// combine EXACTLY into the first quotient, which is the identity
// -y/(ab) + mu/(ac) + (y-mu)/c^2 = (y-mu)^2/(b c^2), so no series is needed
// for them at all.
//
// The two derivations check each other: to leading order the score is
// [y - (y-mu)^2]/(2 th^2) and the Hessian [(y-mu)^2 - y]/th^3, which is its
// derivative.
//
// ⚠️ THE EXPECTED INFORMATION IS NOT FIXED HERE and is measured as wrong: at
// theta = 1e6 it reads -1.7e-16, and an expected information cannot be
// negative.  Its leading order needs one more term of the observed Hessian
// than is written above, the theta^-3 term vanishing under expectation, so it
// is a derivation of its own rather than a transcription of these.

// The three quantities the rewrite needs live in psi_diff.h, shared with
// the beta-binomial, which has the identical shape one family over.

// E[trigamma(Y + theta)] for Y ~ NB2(mu, theta), needed by the expected hessian.
// Sums trigamma(k + theta) * P(Y = k) through the pmf recurrence
// p_{k+1} = p_k * (k + theta) / (k + 1) * mu / (theta + mu), carried in log
// scale until p_k is representable, and stops when the accumulated mass
// reaches 1 - 1e-12 -- the point a far-tail quantile would have located. The
// quantile call an earlier version used here is off limits: this helper runs
// inside d7::par_for workers, and qnbinom's search reaches pbeta, whose
// warning path calls into the R API and killed the process from a worker
// thread on four of the five CI platforms. Everything below is digamma-family
// arithmetic and plain C, which never takes such a path on positive
// arguments. The hard cap covers the geometric tail (decay ratio
// mu/(theta+mu), so ~30(mu+theta)/theta terms reach 1e-12) and the loop
// almost always breaks on the mass long before it.
static double nb_E_trigamma(double mu, double theta) {
    double ratio = mu / (theta + mu);
    double lratio = std::log(mu) - std::log(theta + mu);
    double cap = 100.0 + mu + 20.0 * std::sqrt(mu * (1.0 + mu / theta))
                 + 40.0 * (mu + theta) / theta;
    int kmax = (int) std::min(cap, 2.0e9);

    double s = 0.0, cum = 0.0;
    double lpk = theta * (std::log(theta) - std::log(theta + mu));
    double pk = std::exp(lpk);
    bool logscale = !(pk > 0.0);
    int k = 0;
    for (; k <= kmax; ++k) {
        s += R::trigamma(k + theta) * pk;
        cum += pk;
        if (cum >= 1.0 - 1e-12 && k >= 100) { ++k; break; }
        if (logscale) {
            lpk += lratio + std::log((k + theta) / (k + 1.0));
            pk = std::exp(lpk);
            // leave log scale only once pk is comfortably NORMAL: the first
            // nonzero exp(lpk) is a subnormal with almost no significand, and
            // seeding the multiplicative recurrence there was measured to
            // carry a 2.5x error to the mode (mu = theta = 1e4)
            if (lpk > -640.0) logscale = false;
        } else {
            pk *= (k + theta) / (k + 1.0) * ratio;
        }
    }

    // Tail bound: trigamma is decreasing, so the missing mass contributes at
    // most trigamma(k + theta) * (1 - cum).
    if (cum < 1.0) s += R::trigamma(k + theta) * (1.0 - cum);
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
    double m0 = 0, th0 = 0, th_plus_mu0 = 0;

    if (both_scalar) {
        m0 = mu[0];
        th0 = theta[0];
        th_plus_mu0 = th0 + m0;
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        // LOCAL to the region: a scalar hoisted out of the loop and written
        // inside it is shared once the iterations are split
        double m = m0, th = th0, th_plus_mu = th_plus_mu0;
        if (!both_scalar) {
            m = mu_is_scalar ? mu[0] : mu[i];
            th = theta_is_scalar ? theta[0] : theta[i];
            th_plus_mu = th + m;
        }

        grad_mu[i] = (th / th_plus_mu) * (y[i] / m - 1.0);
        // the two brackets of the rewrite, each cancellation performed
        // symbolically: see the note at the head of this file
        grad_theta[i] = d7::psi_A_rest(y[i], th) +
            d7::psi_Ew((y[i] - m) / th_plus_mu);
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

    double m0 = 0, th0 = 0, th_plus_mu0 = 0, th_plus_mu20 = 0;

    if (both_scalar) {
        m0 = mu[0];
        th0 = theta[0];
        th_plus_mu0 = th0 + m0;
        th_plus_mu20 = th_plus_mu0 * th_plus_mu0;
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        // LOCAL to the region, as above
        double m = m0, th = th0, th_plus_mu = th_plus_mu0,
               th_plus_mu2 = th_plus_mu20;
        if (!both_scalar) {
            m = mu_is_scalar ? mu[0] : mu[i];
            th = theta_is_scalar ? theta[0] : theta[i];
            th_plus_mu = th + m;
            th_plus_mu2 = th_plus_mu * th_plus_mu;
        }

        double res = y[i] - m;

        hess_mu_mu[i] = (y[i] + th) / th_plus_mu2 - y[i] / (m * m);
        // -y/(ab) + mu/(ac) + (y-mu)/c^2 collapses EXACTLY into the first
        // quotient, so those three need no series at all; what is left is
        // the trigamma remainder
        hess_theta_theta[i] = res * res / ((th + y[i]) * th_plus_mu2) +
            d7::psi_T_rest(y[i], th);
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
