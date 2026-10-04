#ifndef D7_PT_NEGBIN2_H
#define D7_PT_NEGBIN2_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "psi_diff.h"

// Negative binomial in the mean and the size, Var = mu + mu^2/theta: one
// function per component and order for the quantities the scalar registry
// reads, called by the vector kernels in negbin.cpp and by the registry in
// d7_ccallable.cpp. negbin.cpp records the rewrite of the theta derivatives
// that keeps them accurate as theta grows; the remainders they read,
//
//   A = psi(y + theta) - psi(theta) - log1p(y/theta)         psi_A_rest
//   E = log1p(w) - w, w = (y - mu)/(theta + mu)              negbin2_Ew
//   T = psi'(y + theta) - psi'(theta) + y/(theta(theta + y)) psi_T_rest
//
// are computed by the caller and passed in.

namespace d7 {

inline double negbin2_score_mu(double y, double m, double th) {
    double th_plus_mu = th + m;
    return (th / th_plus_mu) * (y / m - 1.0);
}

inline double negbin2_score_theta(double A, double E) {
    return A + E;
}

// the score's second bracket, with 1 + w = (y + theta)/(theta + mu) formed
// from its own inputs (psi_diff.h, psi_Ew2)
inline double negbin2_Ew(double y, double m, double th) {
    double th_plus_mu = th + m;
    return d7::psi_Ew2((y + th) / th_plus_mu, (y - m) / th_plus_mu);
}

inline double negbin2_hess_mu_mu(double y, double m, double th) {
    double th_plus_mu = th + m;
    double th_plus_mu2 = th_plus_mu * th_plus_mu;
    return (y + th) / th_plus_mu2 - y / (m * m);
}

// -y/(ab) + mu/(ac) + (y-mu)/c^2 collapses EXACTLY into the first quotient,
// so those three need no series at all; what is left is the trigamma
// remainder T
inline double negbin2_hess_theta_theta(double y, double m, double th,
                                       double T) {
    double th_plus_mu = th + m;
    double th_plus_mu2 = th_plus_mu * th_plus_mu;
    double res = y - m;
    return res * res / ((th + y) * th_plus_mu2) + T;
}

inline double negbin2_expected_mu_mu(double m, double th) {
    double th_plus_mu = th + m;
    return -th / (m * th_plus_mu);
}

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
// The loop is also capped at 1e6 terms. Below theta/mu of about 3e-5 the tail
// needs 1e7 to 1e9 terms, one series then costs 0.2 s or more, and a fit that
// creeps along theta -> 0 pays it at every iteration; that region is the edge
// of the parameter space and its information is not read for an estimate.
inline double negbin2_expected_theta_theta(double mu, double theta) {
    double ratio = mu / (theta + mu);
    double lratio = std::log(mu) - std::log(theta + mu);
    double cap = 100.0 + mu + 20.0 * std::sqrt(mu * (1.0 + mu / theta))
                 + 40.0 * (mu + theta) / theta;
    int kmax = (int) std::min(cap, 1.0e6);
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

// d_mu E_mm = theta(theta + 2 mu)/(mu^2 c^2), with E_mm = -theta/(mu c)
// differentiated in factored form so that nothing cancels as theta -> 0
inline double negbin2_dexpected_mu_mu_mu(double m, double th) {
    double c = th + m, c2 = c * c;
    return th * (th + 2.0 * m) / (m * m * c2);
}

// d_theta of negbin2_expected_theta_theta(). With S = sum_k p_k U_k the
// derivative moves the mass as well as the summand,
//   d_theta S = sum_k p_k (U_theta + U s_theta),
//   U_theta   = -k(2 theta + mu)/(theta^2 c^2) + 2 A3,
//   s_theta   = A1 - log1p(mu/theta) + (mu - k)/c,
// with c = theta + mu and A_r = sum_{j<k} (theta + j)^-r. The mass, its
// seed, its log-scale switch and its stopping rule are those of
// negbin2_expected_theta_theta(); negbin.cpp has the derivation in full.
inline double negbin2_dexpected_theta_theta_theta(double mu, double theta) {
    double ratio = mu / (theta + mu);
    double lratio = std::log(mu) - std::log(theta + mu);
    double cap = 100.0 + mu + 20.0 * std::sqrt(mu * (1.0 + mu / theta))
                 + 40.0 * (mu + theta) / theta;
    int kmax = (int) std::min(cap, 1.0e6);
    const double c = theta + mu, c2 = c * c;
    const double den = theta * c;
    const double L = std::log1p(mu / theta);
    const double th2 = theta * theta;

    double U = 0.0, A1 = 0.0, A3 = 0.0, cum = 0.0;
    double r2 = 0.0;
    double lpk = -theta * std::log1p(mu / theta);
    bool logscale = (lpk <= -640.0);
    double pk = std::exp(lpk);
    for (int k = 0; k <= kmax; ++k) {
        double kd = (double) k;
        double st = A1 - L + (mu - kd) / c;
        double Ut = -kd * (2.0 * theta + mu) / (th2 * c2) + 2.0 * A3;
        r2 += pk * (Ut + U * st);
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
    return r2;
}

inline void negbin2_score_curv(int k, double y, const double* th,
                               double* out) {
    double m = th[0], t = th[1];
    if (k == 0) {
        out[0] = negbin2_score_mu(y, m, t);
        out[1] = negbin2_hess_mu_mu(y, m, t);
    } else {
        out[0] = negbin2_score_theta(d7::psi_A_rest(y, t),
                                     negbin2_Ew(y, m, t));
        out[1] = negbin2_hess_theta_theta(y, m, t, d7::psi_T_rest(y, t));
    }
}

inline void negbin2_info_dinfo(int k, double y, const double* th,
                               double* out) {
    double m = th[0], t = th[1];
    if (k == 0) {
        out[0] = negbin2_expected_mu_mu(m, t);
        out[1] = negbin2_dexpected_mu_mu_mu(m, t);
    } else {
        out[0] = negbin2_expected_theta_theta(m, t);
        out[1] = negbin2_dexpected_theta_theta_theta(m, t);
    }
}

} // namespace d7

#endif
