#include <Rcpp.h>
#include "d7_par.h"
#include "psi_diff.h"
#include "pt_negbin2.h"
#include "pt_sqrt.h"
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

// THE EXPECTED SECOND DERIVATIVE IN theta is negbin2_expected_theta_theta()
// in pt_negbin2.h, with the record of how it is summed.

// THE DERIVATIVES OF E[l_theta_theta] IN (mu, theta).
//
// negbin2_expected_theta_theta() is S(mu, theta) = sum_k p_k U_k with
//   U_k = sum_{j<k} u_j = k/(theta c) - A2_k,   c = theta + mu,
//   A_r,k = sum_{j<k} (theta + j)^-r,
// and differentiating a sum weighted by the mass moves the mass as well as the
// summand, p_k,x = p_k s_x(k) with s the score of the log-mass at k:
//
//   d_x  S = sum p (U_x + U s_x)
//   d_xy S = sum p (U_xy + U_x s_y + U_y s_x + U (s_xy + s_x s_y))
//
// where every piece is an exact polynomial in k or one of the A_r:
//   U_mu      = -k/(theta c^2)
//   U_theta   = -k(2 theta + mu)/(theta^2 c^2) + 2 A3
//   U_mumu    =  2k/(theta c^3)
//   U_mutheta =  k(3 theta + mu)/(theta^2 c^3)
//   U_thth    =  2k(3 theta^2 + 3 theta mu + mu^2)/(theta^3 c^3) - 6 A4
//   s_mu      =  theta (k - mu)/(mu c)
//   s_theta   =  A1 - log1p(mu/theta) + (mu - k)/c
//   s_mumu    =  (k + theta)/c^2 - k/mu^2
//   s_mutheta =  (k - mu)/c^2
//   s_thth    = -A2 + mu/(theta c) + (k - mu)/c^2
//
// The mass comes from the recurrence and the stopping rule of
// negbin2_expected_theta_theta(),
// seed and log-scale switch included, for the reasons recorded there; the tail
// past 1 - 1e-12 of the mass is dropped. nb_dE_ltt1() writes (d_mu S,
// d_theta S) to out in one pass, the d_theta summand being
// negbin2_dexpected_theta_theta_theta_term() of pt_negbin2.h, which the
// registry's d_theta sum also adds; nb_dE_ltt2() writes (d_mumu S,
// d_thth S, d_mutheta S).
static void nb_dE_ltt1(double mu, double theta, double *out) {
    double ratio = mu / (theta + mu);
    double lratio = std::log(mu) - std::log(theta + mu);
    double cap = 100.0 + mu + 20.0 * d7::sqrt_cr(mu * (1.0 + mu / theta))
                 + 40.0 * (mu + theta) / theta;
    int kmax = (int) std::min(cap, 1.0e6);
    const double c = theta + mu, c2 = c * c;
    const double den = theta * c;
    const double L = std::log1p(mu / theta);
    const double th2 = theta * theta;

    double U = 0.0, A1 = 0.0, A3 = 0.0, cum = 0.0;
    double r1 = 0.0, r2 = 0.0;
    double lpk = -theta * std::log1p(mu / theta);
    bool logscale = (lpk <= -640.0);
    double pk = std::exp(lpk);
    for (int k = 0; k <= kmax; ++k) {
        double kd = (double) k;
        double sm = theta * (kd - mu) / (mu * c);
        double Um = -kd / (theta * c2);
        r1 += pk * (Um + U * sm);
        r2 += d7::negbin2_dexpected_theta_theta_theta_term(kd, pk, U, A1, A3, L,
                                                            mu, theta, c, c2, th2);
        cum += pk;
        bool last = (cum >= 1.0 - 1e-12 && k >= 100);
        double tk = theta + kd, iv = 1.0 / tk, iv2 = iv * iv;
        U += (theta * (2.0 * kd - mu) + kd * kd) / (den * tk * tk);
        A1 += iv;
        A3 += iv2 * iv;
        if (last) break;
        if (logscale) {
            lpk += lratio + std::log((k + theta) / (k + 1.0));
            pk = std::exp(lpk);
            if (lpk > -640.0) logscale = false;
        } else {
            pk *= (k + theta) / (k + 1.0) * ratio;
        }
    }
    out[0] = r1;
    out[1] = r2;
}

static void nb_dE_ltt2(double mu, double theta, double *out) {
    double ratio = mu / (theta + mu);
    double lratio = std::log(mu) - std::log(theta + mu);
    double cap = 100.0 + mu + 20.0 * d7::sqrt_cr(mu * (1.0 + mu / theta))
                 + 40.0 * (mu + theta) / theta;
    int kmax = (int) std::min(cap, 1.0e6);
    const double c = theta + mu, c2 = c * c, c3 = c2 * c;
    const double den = theta * c;
    const double L = std::log1p(mu / theta);
    const double th2 = theta * theta, th3 = th2 * theta;

    double U = 0.0, A1 = 0.0, A2 = 0.0, A3 = 0.0, A4 = 0.0, cum = 0.0;
    double r1 = 0.0, r2 = 0.0, r3 = 0.0;
    double lpk = -theta * std::log1p(mu / theta);
    bool logscale = (lpk <= -640.0);
    double pk = std::exp(lpk);
    for (int k = 0; k <= kmax; ++k) {
        double kd = (double) k;
        double sm = theta * (kd - mu) / (mu * c);
        double st = A1 - L + (mu - kd) / c;
        double Um = -kd / (theta * c2);
        double Ut = -kd * (2.0 * theta + mu) / (th2 * c2) + 2.0 * A3;
        double smm = (kd + theta) / c2 - kd / (mu * mu);
        double smt = (kd - mu) / c2;
        double stt = -A2 + mu / (theta * c) + (kd - mu) / c2;
        double Umm = 2.0 * kd / (theta * c3);
        double Umt = kd * (3.0 * theta + mu) / (th2 * c3);
        double Utt = 2.0 * kd * (3.0 * th2 + 3.0 * theta * mu + mu * mu) /
                     (th3 * c3) - 6.0 * A4;
        r1 += pk * (Umm + 2.0 * Um * sm + U * (smm + sm * sm));
        r2 += pk * (Utt + 2.0 * Ut * st + U * (stt + st * st));
        r3 += pk * (Umt + Um * st + Ut * sm + U * (smt + sm * st));
        cum += pk;
        bool last = (cum >= 1.0 - 1e-12 && k >= 100);
        double tk = theta + kd, iv = 1.0 / tk, iv2 = iv * iv;
        U += (theta * (2.0 * kd - mu) + kd * kd) / (den * tk * tk);
        A1 += iv;
        A2 += iv2;
        A3 += iv2 * iv;
        A4 += iv2 * iv2;
        if (last) break;
        if (logscale) {
            lpk += lratio + std::log((k + theta) / (k + 1.0));
            pk = std::exp(lpk);
            if (lpk > -640.0) logscale = false;
        } else {
            pk *= (k + theta) / (k + 1.0) * ratio;
        }
    }
    out[0] = r1;
    out[1] = r2;
    out[2] = r3;
}

// The derivatives of the expected information in the parameters.
//   E_mm = -theta/(mu c) = -1/mu + 1/c,  E_mt = 0,  E_tt = negbin2_expected_theta_theta(),
// with E_mm's written in factored form so that nothing cancels as theta -> 0:
//   d_mu E_mm      =  theta(theta + 2 mu)/(mu^2 c^2),  d_theta E_mm = -1/c^2
//   d_mumu E_mm    = -2 theta (c^2 + c mu + mu^2)/(mu^3 c^3)
//   d_mutheta E_mm =  2/c^3,  d_thth E_mm = 2/c^3.
// [[Rcpp::export]]
List negbin_dexpected1_cpp(NumericVector y, NumericVector mu, NumericVector theta, int threads = 1) {
    int n = y.size();
    bool m_s = (mu.size() == 1), t_s = (theta.size() == 1);
    const double *mp = mu.begin(), *tp = theta.begin();
    NumericVector zero(n);
    NumericVector mm_m(n), mm_t(n), tt_m(n), tt_t(n);
    double *a = mm_m.begin(), *b = mm_t.begin(), *e = tt_m.begin(),
           *f = tt_t.begin();
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m_s ? mp[0] : mp[i];
        double th = t_s ? tp[0] : tp[i];
        double c = th + m, c2 = c * c;
        a[i] = d7::negbin2_dexpected_mu_mu_mu(m, th);
        b[i] = -1.0 / c2;
        double r[2];
        nb_dE_ltt1(m, th, r);
        e[i] = r[0];
        f[i] = r[1];
    });
    return List::create(
        Named("mu_mu_mu") = mm_m, Named("mu_mu_theta") = mm_t,
        Named("theta_theta_mu") = tt_m, Named("theta_theta_theta") = tt_t,
        Named("mu_theta_mu") = zero, Named("mu_theta_theta") = clone(zero));
}

// [[Rcpp::export]]
List negbin_dexpected2_cpp(NumericVector y, NumericVector mu, NumericVector theta, int threads = 1) {
    int n = y.size();
    bool m_s = (mu.size() == 1), t_s = (theta.size() == 1);
    const double *mp = mu.begin(), *tp = theta.begin();
    NumericVector zero(n);
    NumericVector mm_mm(n), mm_tt(n), mm_mt(n), tt_mm(n), tt_tt(n), tt_mt(n);
    double *a = mm_mm.begin(), *b = mm_tt.begin(), *cc = mm_mt.begin(),
           *e = tt_mm.begin(), *f = tt_tt.begin(), *g = tt_mt.begin();
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m_s ? mp[0] : mp[i];
        double th = t_s ? tp[0] : tp[i];
        double c = th + m, c2 = c * c, c3 = c2 * c;
        double r[3];
        nb_dE_ltt2(m, th, r);
        a[i] = -2.0 * th * (c2 + c * m + m * m) / (m * m * m * c3);
        b[i] = 2.0 / c3;
        cc[i] = 2.0 / c3;
        e[i] = r[0];
        f[i] = r[1];
        g[i] = r[2];
    });
    return List::create(
        Named("mu_mu_mu_mu") = mm_mm, Named("mu_mu_theta_theta") = mm_tt,
        Named("mu_mu_mu_theta") = mm_mt,
        Named("theta_theta_mu_mu") = tt_mm, Named("theta_theta_theta_theta") = tt_tt,
        Named("theta_theta_mu_theta") = tt_mt,
        Named("mu_theta_mu_mu") = zero, Named("mu_theta_theta_theta") = clone(zero),
        Named("mu_theta_mu_theta") = clone(zero));
}

// [[Rcpp::export]]
List negbin_gradient_cpp(NumericVector y, NumericVector mu, NumericVector theta,
                        int threads = 1) {
    int n = y.size();
    NumericVector grad_mu(n);
    NumericVector grad_theta(n);

    bool mu_is_scalar = (mu.size() == 1);
    bool theta_is_scalar = (theta.size() == 1);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double th = theta_is_scalar ? theta[0] : theta[i];

        grad_mu[i] = d7::negbin2_score_mu(y[i], m, th);
        // the two brackets of the rewrite, each cancellation performed
        // symbolically: see the note at the head of this file
        grad_theta[i] = d7::negbin2_score_theta(d7::psi_A_rest(y[i], th),
                                                d7::negbin2_Ew(y[i], m, th));
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
        // LOCAL to the region: a scalar hoisted out of the loop and written
        // inside it is shared once the iterations are split
        double m = m0, th = th0, th_plus_mu = th_plus_mu0,
               th_plus_mu2 = th_plus_mu20;
        if (!both_scalar) {
            m = mu_is_scalar ? mu[0] : mu[i];
            th = theta_is_scalar ? theta[0] : theta[i];
            th_plus_mu = th + m;
            th_plus_mu2 = th_plus_mu * th_plus_mu;
        }

        double res = y[i] - m;

        hess_mu_mu[i] = d7::negbin2_hess_mu_mu(y[i], m, th);
        hess_theta_theta[i] = d7::negbin2_hess_theta_theta(
            y[i], m, th, d7::psi_T_rest(y[i], th));
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
        hmm0 = d7::negbin2_expected_mu_mu(m, th);
        htt0 = d7::negbin2_expected_theta_theta(m, th);
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double hmm = hmm0, htt = htt0;
        if (!both_scalar) {
            double m = mu_is_scalar ? mu[0] : mu[i];
            double th = theta_is_scalar ? theta[0] : theta[i];
            hmm = d7::negbin2_expected_mu_mu(m, th);
            htt = d7::negbin2_expected_theta_theta(m, th);
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
