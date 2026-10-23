#ifndef D7_PT_NEGBIN2_H
#define D7_PT_NEGBIN2_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "psi_diff.h"
#include "pt_sqrt.h"
#include "pt_discrete_pow.h"

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

// THE DERIVATIVES IN theta AT y = 0 AND y = 1. There the general forms above
// lose relative accuracy as theta^2/mu: at y = 1 the score is of order
// mu/theta^2 and is formed as A + E from two terms of order 1/theta (1.4e-9
// at mu = 1e-6, theta = 5; 4.8e-6 at mu = 1e-8, theta = 50, and 2.8e-5 in the
// Hessian), and at y = 0 orders three and four lose the same way. With
// h = theta L and L = log1p(mu/theta), c = theta + mu,
//   log f(0) = -h,   log f(1) = log mu - h - L,
// and only h' = L - x/(1 + x), x = mu/theta, cancels: below x = 0.5 it is the
// series sum_{n>=2} (-1)^n (n - 1)/n x^n. The higher derivatives of h and
// every derivative of L are sums of terms of one sign,
//   h''   = -mu^2/(theta c^2),
//   h'''  = mu^2/(theta^2 c^2) + 2 mu^2/(theta c^3),
//   h'''' = -2 mu^2/(theta^3 c^2) - 4 mu^2/(theta^2 c^3) - 6 mu^2/(theta c^4),
//   L' = -mu/(theta c),  L'' = mu (2 theta + mu)/(theta^2 c^2),
//   L''' = -2 mu (3 theta^2 + 3 theta mu + mu^2)/(theta^3 c^3),
//   L'''' = 6 mu (4 theta^3 + 6 theta^2 mu + 4 theta mu^2 + mu^3)/(theta^4 c^4),
// written below in 1/theta and 1/c so that nothing overflows at large theta.
// Against mpmath (misure_gas_congiunto/c/fase1/c/nb2_y01_proto.py) at
// 160 digits, mu 1e-8 to 1e3 and theta 1e-2 to 1e12: 1e-12 at worst (mu =
// 1e3, theta = 2, order three), 4e-16 on the small-mu grid of item (c).
inline double negbin2_h1(double x) {
    if (x < 0.5) {
        double s = 0.0, p = x * x;
        for (int n = 2; n < 200; ++n) {
            const double t = p * (n - 1.0) / n;
            s += (n % 2 == 0) ? t : -t;
            if (t < 1e-17 * x * x * 0.5) break;
            p *= x;
        }
        return s;
    }
    return std::log1p(x) - x / (1.0 + x);
}

// d^order/dtheta^order log f(y) for y = 0 or 1, order 1 to 4
inline double negbin2_theta_y01(double y, double m, double th, int order) {
    // with u = 1/theta and v = 1/c the derivatives of L are, from
    // L^(k) = (-1)^k (k - 1)! (u^k - v^k) and u - v = a v,
    //   L' = -a v, L'' = a v (u + v), L''' = -2 a v (u^2 + u v + v^2),
    //   L'''' = 6 a v (u^3 + u^2 v + u v^2 + v^3),
    // and those of h are a b v times a polynomial in u and v
    const double c = th + m, a = m / th, b = m / c;
    const double u = 1.0 / th, v = 1.0 / c, av = a * v, abv = a * b * v;
    double h, L;
    switch (order) {
    case 1:
        h = negbin2_h1(a);
        L = -av;
        break;
    case 2:
        h = -abv;
        L = av * (u + v);
        break;
    case 3:
        h = abv * (u + 2.0 * v);
        L = -2.0 * av * (u * u + u * v + v * v);
        break;
    default:
        h = -abv * (2.0 * u * u + 4.0 * u * v + 6.0 * v * v);
        L = 6.0 * av * (u * u * u + u * u * v + u * v * v + v * v * v);
    }
    return (y < 0.5) ? -h : -h - L;
}

// the score in theta at one observation, through the closed forms at y = 0
// and y = 1 and the two remainders elsewhere
inline double negbin2_score_theta_y(double y, double m, double th) {
    if (y < 1.5) return negbin2_theta_y01(y, m, th, 1);
    return negbin2_score_theta(d7::psi_A_rest(y, th), negbin2_Ew(y, m, th));
}

inline double negbin2_hess_mu_mu(double y, double m, double th) {
    // -T_2 of count_bracket() (pt_discrete_pow.h): written directly, the two
    // quotients agree to a relative theta/mu and the value, of order theta,
    // loses that many digits at a small size
    return -count_bracket(y, m, th, 2);
}

// -y/(ab) + mu/(ac) + (y-mu)/c^2 collapses EXACTLY into the first quotient,
// so those three need no series at all; what is left is the trigamma
// remainder T
inline double negbin2_hess_theta_theta(double y, double m, double th,
                                       double T) {
    if (y < 1.5) return negbin2_theta_y01(y, m, th, 2);
    double th_plus_mu = th + m;
    double th_plus_mu2 = th_plus_mu * th_plus_mu;
    double res = y - m;
    return res * res / ((th + y) * th_plus_mu2) + T;
}

inline double negbin2_expected_mu_mu(double m, double th) {
    double th_plus_mu = th + m;
    return -th / (m * th_plus_mu);
}

// THE STOPPING RULE OF THE SERIES BELOW, shared by every sum over the support
// of the negative binomial's expected information and its derivatives
// (negbin.cpp, negbin_hd.cpp). The loop stops past the mode, and not before
// k = 100, once p_k (1 + k)^2 falls below 1e-17, (1 + k)^2 bounding the
// growth of the largest summand. An earlier rule stopped when the accumulated
// mass reached 1 - 1e-12 and added U (1 - cum) for the tail; measured against
// mpmath on 31 cells with theta < 1e4, that left 5e-11 on E[l_theta_theta],
// 1.5e-7 on its derivative in mu and 3.7e-9 on its derivative in theta, and
// this rule leaves 2.6e-14, 6.3e-13 and 1.9e-11 for 1.7 times the terms. The
// accumulated mass cannot serve at a tighter tolerance: its own rounding keeps
// it below 1 - 1e-16. The cap allows the geometric tail, of ratio
// mu/(theta + mu), the extra terms down to 1e-17.
inline bool negbin2_series_done(double kd, double pk, double mode) {
    return kd >= 100.0 && kd > mode && pk * (1.0 + kd) * (1.0 + kd) <= 1e-17;
}

inline double negbin2_series_mode(double mu, double theta) {
    return (theta > 1.0) ? mu * (theta - 1.0) / theta : 0.0;
}

inline int negbin2_series_cap(double mu, double theta) {
    const double cap = 100.0 + mu + 20.0 * d7::sqrt_cr(mu * (1.0 + mu / theta))
                       + 80.0 * (mu + theta) / theta;
    return (int) std::min(cap, 1.0e6);
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
// the loop stops by negbin2_series_done().  The quantile call an earlier version
// used here is off limits: this helper runs inside d7::par_for workers, and
// qnbinom's search reaches pbeta, whose warning path calls into the R API and
// killed the process from a worker thread on four of the five CI platforms.
// Everything below is plain C arithmetic, which never takes such a path.  The
// hard cap covers the geometric tail (negbin2_series_cap()) and the loop
// almost always stops by its rule long before it.
// The loop is also capped at 1e6 terms. Below theta/mu of about 3e-5 the tail
// needs 1e7 to 1e9 terms, one series then costs 0.2 s or more, and a fit that
// creeps along theta -> 0 pays it at every iteration; that region is the edge
// of the parameter space and its information is not read for an estimate.
inline double negbin2_expected_theta_theta(double mu, double theta) {
    double ratio = mu / (theta + mu);
    double lratio = std::log(mu) - std::log(theta + mu);
    const int kmax = negbin2_series_cap(mu, theta);
    const double mode = negbin2_series_mode(mu, theta);
    const double den = theta * (theta + mu);

    double s = 0.0, U = 0.0;
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
        double kd = (double) k, tk = theta + kd;
        bool last = negbin2_series_done(kd, pk, mode);
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
//
// THE RUNNING FORMS OF s_theta, s_thetatheta, U_theta AND U_thetatheta.
// Written from the closed forms above, U_theta = -k(2 theta + mu)/(theta^2 c^2)
// + 2 A3 and U_thetatheta = 2k(3 theta^2 + 3 theta mu + mu^2)/(theta^3 c^3)
// - 6 A4 are each a difference of two sums of the same size, and so are
//   s_theta      = A1 - log1p(mu/theta) + (mu - k)/c,
//   s_thetatheta = -A2 + mu/(theta c) + (k - mu)/c^2.
// Their loss grows with theta, because the information tends to zero with
// theta while the pieces do not. Each is
// here the running sum of its merged increment, which is the same difference
// taken term by term (each A_r absorbs its constant, as the merged summand of
// U does and as nb_E_dtheta of negbin_hd.cpp does for U_theta and
// U_thetatheta):
//   s_theta,k+1      - s_theta,k      = (mu - k)/(t c)
//   s_thetatheta,k+1 - s_thetatheta,k = (k - mu)(t + c)/(t^2 c^2)
//   U_theta,k+1      - U_theta,k      = [theta^2 mu (3 theta + 2 mu)
//                         - (2 theta + mu) k (3 theta^2 + 3 k theta + k^2)]
//                         / (t^3 theta^2 c^2)
//   U_thetatheta,k+1 - U_thetatheta,k = [-2 theta^3 mu (6 theta^2 + 8 theta mu
//                         + 3 mu^2) + 2 K k (4 theta^3 + 6 k theta^2
//                         + 4 k^2 theta + k^3)] / (t^4 theta^3 c^3),
// t = theta + k and K = 3 theta^2 + 3 theta mu + mu^2. All four vanish at
// k = 0 except the two s, whose values there are
//   s_theta,0      = mu/c - log1p(mu/theta)  (the series negbin2_h1() below
//                    mu/theta = 0.5, where the two terms cancel),
//   s_thetatheta,0 = mu^2/(theta c^2).
// The loop starts at k = 0, so no remainder is needed to seed U_theta and
// U_thetatheta. Every increment is formed from iv = 1/t by multiplication, so
// that the four sums cost no division beyond the one iv the loop already took.
// Against nb_ref.json (misure_gas_congiunto/c/fase3/N2), over the negbin2
// cells with mu from 1e-6 to 1e3 and theta from 0.05 to 1e5: the worst error
// of the first derivatives falls from 4.4e-9 to 2.6e-11 and that of the
// second derivatives from 3.4e-8 to 9.8e-10; what remains of the second is
// d_mumu E_tt, which reads neither U_theta nor s_theta.
struct Negbin2ThetaRun {
    double st, stt, Ut, Utt;   // the four quantities at the current k
    double mu, theta, c, invc, ic2, ic3, ic4, th2, th3, A3, B3, A4, K4;
    // ratio = mu/(theta + mu) and L = log1p(mu/theta), formed by the caller
    Negbin2ThetaRun(double mu_, double theta_, double ratio, double L)
        : Ut(0.0), Utt(0.0), mu(mu_), theta(theta_) {
        c = theta + mu;
        invc = 1.0 / c;
        const double c2 = c * c, c3 = c2 * c;
        ic2 = invc * invc;
        th2 = theta * theta;
        th3 = th2 * theta;
        ic3 = 1.0 / (th2 * c2);
        ic4 = 1.0 / (th3 * c3);
        A3 = th2 * mu * (3.0 * theta + 2.0 * mu);
        B3 = 2.0 * theta + mu;
        A4 = -2.0 * th3 * mu * (6.0 * th2 + 8.0 * theta * mu + 3.0 * mu * mu);
        K4 = 3.0 * th2 + 3.0 * theta * mu + mu * mu;
        const double x = mu / theta;
        st = (x < 0.5) ? -negbin2_h1(x) : ratio - L;
        stt = mu * mu / (theta * c2);
    }
    // the increments of s_theta and U_theta at k = kd, tk = theta + kd, iv = 1/tk
    inline void step1(double kd, double iv) {
        const double iv2 = iv * iv;
        st += (mu - kd) * iv * invc;
        Ut += (A3 - B3 * kd * (3.0 * th2 + 3.0 * kd * theta + kd * kd))
              * (iv2 * iv) * ic3;
    }
    // the increments of s_thetatheta and U_thetatheta
    inline void step2(double kd, double tk, double iv) {
        const double iv2 = iv * iv;
        stt += (kd - mu) * (tk + c) * iv2 * ic2;
        Utt += (A4 + 2.0 * K4 * kd
                * (4.0 * th3 + 6.0 * kd * th2 + 4.0 * kd * kd * theta
                   + kd * kd * kd)) * (iv2 * iv2) * ic4;
    }
};

// The summand at k, from U and the running s_theta and U_theta at k, is its
// own function: the vector kernel, which forms d_mu S in the same pass, adds
// it too.
inline double negbin2_dexpected_theta_theta_theta_term(double pk, double U,
                                                       double st, double Ut) {
    return pk * (Ut + U * st);
}

inline double negbin2_dexpected_theta_theta_theta(double mu, double theta) {
    double ratio = mu / (theta + mu);
    double lratio = std::log(mu) - std::log(theta + mu);
    const int kmax = negbin2_series_cap(mu, theta);
    const double mode = negbin2_series_mode(mu, theta);
    const double den = theta * (theta + mu);
    const double L = std::log1p(mu / theta);
    Negbin2ThetaRun run(mu, theta, ratio, L);

    double U = 0.0;
    double r2 = 0.0;
    double lpk = -theta * std::log1p(mu / theta);
    bool logscale = (lpk <= -640.0);
    double pk = std::exp(lpk);
    for (int k = 0; k <= kmax; ++k) {
        double kd = (double) k;
        r2 += negbin2_dexpected_theta_theta_theta_term(pk, U, run.st, run.Ut);
        bool last = negbin2_series_done(kd, pk, mode);
        double tk = theta + kd, iv = 1.0 / tk;
        U += (theta * (2.0 * kd - mu) + kd * kd) / (den * tk * tk);
        run.step1(kd, iv);
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
        out[0] = negbin2_score_theta_y(y, m, t);
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

// The mixed observed entry, l_mu_theta = (y - mu)/(theta + mu)^2, as
// negbin_hessian_cpp() forms it.
inline double negbin2_hess_mu_theta(double y, double m, double th) {
    double th_plus_mu = th + m;
    double th_plus_mu2 = th_plus_mu * th_plus_mu;
    double res = y - m;
    return res / th_plus_mu2;
}

// d_theta E_mm = -1/c^2 with c = theta + mu, as negbin_dexpected1_cpp() forms it.
inline double negbin2_dexpected_mu_mu_theta(double m, double th) {
    double c = th + m, c2 = c * c;
    return -1.0 / c2;
}

// d_mu E_tt, the mu-derivative of the expected dispersion entry. It is the
// first component of nb_dE_ltt1() in negbin.cpp, with the same mass
// recurrence, seed, log-scale switch and stopping rule as
// negbin2_expected_theta_theta(). The summand is
// p_k (U_mu + U s_mu) with U_mu = -k/(theta c^2), s_mu = theta (k - mu)/(mu c).
inline double negbin2_dexpected_theta_theta_mu(double mu, double theta) {
    double ratio = mu / (theta + mu);
    double lratio = std::log(mu) - std::log(theta + mu);
    const int kmax = negbin2_series_cap(mu, theta);
    const double mode = negbin2_series_mode(mu, theta);
    const double c = theta + mu, c2 = c * c;
    const double den = theta * c;

    double U = 0.0;
    double r1 = 0.0;
    double lpk = -theta * std::log1p(mu / theta);
    bool logscale = (lpk <= -640.0);
    double pk = std::exp(lpk);
    for (int k = 0; k <= kmax; ++k) {
        double kd = (double) k;
        double sm = theta * (kd - mu) / (mu * c);
        double Um = -kd / (theta * c2);
        r1 += pk * (Um + U * sm);
        bool last = negbin2_series_done(kd, pk, mode);
        double tk = theta + kd;
        U += (theta * (2.0 * kd - mu) + kd * kd) / (den * tk * tk);
        if (last) break;
        if (logscale) {
            lpk += lratio + std::log((k + theta) / (k + 1.0));
            pk = std::exp(lpk);
            if (lpk > -640.0) logscale = false;
        } else {
            pk *= (k + theta) / (k + 1.0) * ratio;
        }
    }
    return r1;
}

// The routers of the joint registry, one per order, for (mu, theta) in that
// order. Each fills a full P x P matrix (P = 2) column-major, out[a + P b],
// and the derivative of the expected information a P x P x P array, out[a +
// P b + P P c] being the derivative of E[l_ab] in parameter c.
inline void negbin2_joint_grad(double y, const double* th, double* out) {
    double m = th[0], t = th[1];
    out[0] = negbin2_score_mu(y, m, t);
    out[1] = negbin2_score_theta_y(y, m, t);
}

inline void negbin2_joint_hess(double y, const double* th, double* out) {
    double m = th[0], t = th[1];
    out[0] = negbin2_hess_mu_mu(y, m, t);
    out[3] = negbin2_hess_theta_theta(y, m, t, d7::psi_T_rest(y, t));
    out[1] = out[2] = negbin2_hess_mu_theta(y, m, t);
}

inline void negbin2_joint_info(double y, const double* th, double* out) {
    double m = th[0], t = th[1];
    out[0] = negbin2_expected_mu_mu(m, t);
    out[3] = negbin2_expected_theta_theta(m, t);
    out[1] = out[2] = 0.0;
}

inline void negbin2_joint_dinfo(double y, const double* th, double* out) {
    double m = th[0], t = th[1];
    for (int j = 0; j < 8; ++j) out[j] = 0.0;
    out[0] = negbin2_dexpected_mu_mu_mu(m, t);
    out[4] = negbin2_dexpected_mu_mu_theta(m, t);
    out[3] = negbin2_dexpected_theta_theta_mu(m, t);
    out[7] = negbin2_dexpected_theta_theta_theta(m, t);
}

} // namespace d7

#endif
