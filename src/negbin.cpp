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
// The three quantities the rewrite needs live in psi_diff.h, shared with
// the beta-binomial, which has the identical shape one family over.

// THE EXPECTED SECOND DERIVATIVE IN theta, E[l_theta_theta], AS ONE SUM.
//
// It used to be assembled from two pieces -- E[psi'(Y+theta)] with psi'(theta)
// subtracted at the call site, plus mu/(theta(theta+mu)) -- and each of those
// carried a cancellation of its own.  Both are removed here.
//
// FIRST, the polygamma difference.  The shift is Y, which is a count, and for
// an INTEGER shift the difference has an exact recurrence whose terms carry
// one sign,
//
//   psi'(x + k) - psi'(x) = - sum_{j<k} 1/(x + j)^2,
//
// so nothing of the leading behaviour is formed and then subtracted, and the
// term costs one division where it used to cost a call to trigamma.  Measured
// over 277 cells spanning mu from 0.1 to 1e5 and theta from 0.05 to 1e6, the
// spelling this replaces is 2.08e-03 out at mu = 0.1, theta = 5.012e5 against
// a sum taken independently in R; on the 57600 rows of a real fit it cost
// 3.570 s an evaluation against 0.070 s here.  The same identity serves every
// order and negbin_hd.cpp uses it at two more.
//
// SECOND, the two pieces themselves.  They are each of order mu/theta^2 while
// their sum is of order mu^2/(2 theta^4), so composing them loses thirteen
// digits at theta = 1e6 and the result reads NEGATIVE past it -- measured at
// mu = 100, -4.47e-21 at theta = 1e7 -- which an expected information cannot
// be.  The second piece is itself an expectation over the same support, since
// E[sum_{j<Y} c] = c E[Y] = c mu with c = 1/(theta(theta+mu)), so the two merge
// into ONE accumulation whose term is
//
//   c - 1/(theta+j)^2 = (theta(2j - mu) + j^2) / (theta (theta+mu) (theta+j)^2),
//
// a quotient of exact polynomials, and it costs nothing.  ⚠️ That is measured
// BACK TO BACK, which is the only way it can be: two whole fits of the same
// model run one after the other, the composition at 263.6 s and this form at
// 264.8 s with the outer trajectory identical evaluation for evaluation.  Runs
// taken at DIFFERENT moments of the same session read 270.3 and 271.3 against
// 292.8 and 292.2, which looks like eight per cent and is the machine: per
// term the two forms measure within 1.4% of each other at every mean from 1 to
// 5e4, with this one never the slower.
//
// ⚠️ That does not remove the cancellation, it moves it from order theta to
// order theta/mu: the terms are of size mu/theta^3 and change sign at
// j = mu/2, while the sum is of size mu^2/(2 theta^4).  What it buys is
// measured, against an asymptote DERIVED rather than fitted -- expanding both
// pieces in 1/theta, the theta^-3 term of the sum contributes mu^2/theta^4 and
// the theta^-4 term -3mu^2/(2 theta^4), so the information tends to
// mu^2/(2 theta^4).  Over mu in {1, 4, 100} and theta from 10 to 1e8 this form
// is POSITIVE at every one of 24 cells and its distance from that asymptote
// FALLS monotonically (at mu = 1: 5.22e-06, 3.02e-07, 1.79e-08 at theta = 1e6,
// 1e7, 1e8), while the composition it replaces reads -1.63e-25 at mu = 4,
// theta = 1e7 and -5.63e-28 at mu = 100, theta = 1e8.
//
// ⚠️ Past theta = 1e5 a reference summed in R DIVERGES from this form -- its
// own distance from the asymptote grows, 2.02e-05 to 1.07e-01 over the same
// three points -- so it is the reference that fails there and not the kernel.
// And the Poisson limit is NOT a reference for it at all: under a Poisson mass
// the answer is exactly 3 mu^2/(2 theta^4), a factor of three, because
// E[Y(Y-1)] is mu^2 there against mu^2(1 + 1/theta) here and the difference is
// precisely the mu^2/theta^4 term the derivation above turns on.  That factor
// coming out as exactly 3 is what confirms the derivation.
//
// The mass comes from the pmf recurrence p_{k+1} = p_k * (k + theta) / (k + 1)
// * mu / (theta + mu), carried in log scale until p_k is representable, and
// the loop stops when the accumulated mass reaches 1 - 1e-12 -- the point a
// far-tail quantile would have located.  The quantile call an earlier version
// used here is off limits: this helper runs inside d7::par_for workers, and
// qnbinom's search reaches pbeta, whose warning path calls into the R API and
// killed the process from a worker thread on four of the five CI platforms.
// Everything below is plain C arithmetic, which never takes such a path.  The
// hard cap covers the geometric tail (decay ratio mu/(theta+mu), so
// ~30(mu+theta)/theta terms reach 1e-12) and the loop almost always breaks on
// the mass long before it.
static double nb_E_ltt(double mu, double theta) {
    double ratio = mu / (theta + mu);
    double lratio = std::log(mu) - std::log(theta + mu);
    double cap = 100.0 + mu + 20.0 * std::sqrt(mu * (1.0 + mu / theta))
                 + 40.0 * (mu + theta) / theta;
    int kmax = (int) std::min(cap, 2.0e9);
    const double den = theta * (theta + mu);

    double s = 0.0, cum = 0.0, U = 0.0;
    // log P(Y = 0) is -theta log1p(mu/theta) and NOT theta (log theta -
    // log(theta + mu)): at theta = 1.585e5 with mu = 0.1 those two logarithms
    // are 11.9736 apiece and their difference 6.31e-07, so the seed loses
    // seven digits and the mass recurrence carries the loss to every term.
    // Measured with nothing else changed, that spelling reads 3.52e-07 out
    // where log1p reads 5.74e-13.
    double lpk = -theta * std::log1p(mu / theta);
    // The log scale is left on the SAME threshold the loop switches at, and
    // not on exp(lpk) merely being nonzero. Between the smallest denormal and
    // the smallest normal double the exponential returns a number with almost
    // no significand, so a rule reading the value leaves log scale exactly
    // where the multiplicative recurrence must not be seeded. Measured at
    // mu = 1.702e5, theta = 100, where log P(Y = 0) is -744.0 and the mass at
    // zero is 5e-324: this expectation came back 3.7e-02 out, and the error
    // runs smoothly up to that edge (2.3e-06 at -734, 1.4e-03 at -739).
    bool logscale = (lpk <= -640.0);
    double pk = std::exp(lpk);
    int k = 0;
    for (; k <= kmax; ++k) {
        s += U * pk;
        cum += pk;
        bool last = (cum >= 1.0 - 1e-12 && k >= 100);
        double kd = (double) k, tk = theta + kd;
        U += (theta * (2.0 * kd - mu) + kd * kd) / (den * tk * tk);
        if (last) { ++k; break; }
        if (logscale) {
            lpk += lratio + std::log((k + theta) / (k + 1.0));
            pk = std::exp(lpk);
            if (lpk > -640.0) logscale = false;
        } else {
            pk *= (k + theta) / (k + 1.0) * ratio;
        }
    }
    if (cum < 1.0) s += U * (1.0 - cum);
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
            d7::psi_Ew2((y[i] + th) / th_plus_mu, (y[i] - m) / th_plus_mu);
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
        htt0 = nb_E_ltt(m, th);
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double hmm = hmm0, htt = htt0;
        if (!both_scalar) {
            double m = mu_is_scalar ? mu[0] : mu[i];
            double th = theta_is_scalar ? theta[0] : theta[i];
            double th_plus_mu = th + m;
            hmm = -th / (m * th_plus_mu);
            htt = nb_E_ltt(m, th);
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
