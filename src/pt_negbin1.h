#ifndef D7_PT_NEGBIN1_H
#define D7_PT_NEGBIN1_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "psi_diff.h"
#include "pt_sqrt.h"

// The negative binomial with a variance linear in the mean, in (mu, theta),
// r = mu/theta: the parts of the log-mass (nb1_parts), the support sums that
// carry its expected information (nb1_E_Pr, nb1_G_derivs1), and one
// function per component and order for the quantities the scalar registry
// reads, called by the vector kernels in negbin1.cpp, which records the
// derivation, and by the registry (d7_ccallable.cpp).

struct NB1parts {
    double P, Q, Pr, Pth, Qth;
    double rm, rt, rmt, rtt;      // derivatives of r in (mu, theta)
};

inline NB1parts nb1_parts(double y, double mu, double th, int order) {
    NB1parts z;
    double r = mu / th, om = 1.0 + th;
    // P and Pr both vanish at the Poisson limit theta -> 0, where r runs
    // away, and the consumers below divide them by theta^2 and theta^4.
    // Written directly the score is 1.5 per cent out at theta = 1e-6, a
    // factor of 400 out at 1e-8 and of the WRONG SIGN at 1e-10.
    //
    // For P the two logarithms combine exactly: with y/r = y theta/mu,
    //   log1p(y/r) - log1p(theta) = log1p( theta (y - mu) / (mu (1 + th)) ),
    // so nothing of the leading behaviour is formed and then subtracted.
    z.P   = d7::psi_A_rest(y, r) + std::log1p(th * (y - mu) / (mu * om));
    z.Q   = -r / om + y / th - y / om;
    z.rm  = 1.0 / th;
    z.rt  = -mu / (th * th);
    if (order < 2) return z;
    // Pr is the only place psi' enters, and below psi_T_rest's crossover that
    // is two trigamma calls an observation. The gradient does not read it, so
    // the second-order block is written only when it is asked for.
    z.Pr  = d7::psi_T_rest(y, r) - (y / r) / (y + r);
    z.Pth = -1.0 / om;
    z.Qth = r / (om * om) - y / (th * th) + y / (om * om);
    z.rmt = -1.0 / (th * th);
    z.rtt = 2.0 * mu / (th * th * th);
    return z;
}

// E[psi'(Y + r) - psi'(r)] under the family itself, needed by the expected
// hessian. The same device the NB2 kernel uses: there is no closed form, and
// a series against the exact mass is better than a quadrature.
//
// It is the DIFFERENCE that is summed, term by term, and not the expectation
// of psi'(Y + r) with psi'(r) taken off at the end: the two agree to leading
// order as r runs away, so subtracting after summing loses the answer where
// subtracting inside the sum costs nothing.
//
// The term is that difference at an INTEGER shift, so it has an exact
// recurrence whose terms carry one sign,
//
//   psi'(r + k) - psi'(r) = - sum_{j<k} 1/(r + j)^2,
//
// and it is accumulated beside the mass -- T_0 = 0, T_{k+1} = T_k -
// 1/(r + k)^2 -- rather than evaluated by psi_T_rest at every term. Below
// r = 100 that helper takes two trigamma calls per term, the second always at
// the same argument: measured with the repetition loop sized by elapsed time,
// 1451 ns a term at r = 0.5 and 472 ns at r = 6 against 17.4 ns and 4.5 ns
// here, and above r = 100, where psi_T_rest is already its asymptotic series,
// 8.2 ns against 2.2 ns. On the 57600 fitted means of a real design one
// evaluation goes from 0.360 s to 0.060 s at theta = 0.2 and from 40.14 s to
// 0.25 s at theta = 100, agreeing with the form it replaces to 1.0e-13.
//
// The mass comes from the recurrence p_{k+1} = p_k (k + r)/(k + 1) *
// theta/(1 + theta), carried in log scale until it is representable, and the
// loop stops once the accumulated mass reaches 1 - 1e-12. That is the point
// the far-tail quantile an earlier version called would have located, and
// the call was R::qnbinom, which is off limits here for the reason
// nb_E_trigamma states one family over: this helper runs inside d7::par_for
// workers, and qnbinom's search reaches pbeta, whose warning path calls into
// the R API and killed the process from a worker thread on four of five CI
// platforms. Removing it is what lets the caller run in parallel at all, and
// it was also the cost. Measured on a 28800-cell design with an offset,
// where every observation carries its own mean: one evaluation of the
// expected information cost 1.3133 s and costs 0.3456 s summed this way,
// 0.0468 s over eight threads, the two routes agreeing to 1.5e-11. The old
// form summed at least a hundred terms whatever the mass did, so on a design
// of small counts the gain is larger.
//
// The hard cap covers the geometric tail, whose decay ratio is
// theta/(1 + theta); the loop almost always breaks on the mass long before.
inline double nb1_E_Pr(double mu, double th) {
    double r = mu / th;
    double ratio = th / (1.0 + th);
    double lratio = std::log(th) - std::log1p(th);
    double cap = 100.0 + mu + 20.0 * d7::sqrt_cr(mu * (1.0 + th))
                 + 40.0 / (-std::log(ratio));
    int kmax = (int) std::min(cap, 2.0e9);

    double lpk = -r * std::log1p(th);      // log P(Y = 0), which is r log(1 - ratio)
    // The log scale is left on the SAME threshold the loop switches at, and
    // not on exp(lpk) merely being nonzero: at mu = 918.832, theta = 0.5 the
    // size is 1837.7 and lpk is -745.1, so P(Y = 0) is a subnormal with
    // almost no significand. Seeding the multiplicative recurrence there put
    // mu_theta at -2.73e-02 where it is 1.2089e-04.
    bool logscale = (lpk <= -640.0);
    double pk = std::exp(lpk);
    double s = 0.0, cum = 0.0, T = 0.0;
    for (int k = 0; k <= kmax; ++k) {
        double kd = (double) k;
        s += T * pk;
        cum += pk;
        if (cum >= 1.0 - 1e-12) break;
        double v = 1.0 / (r + kd);
        T -= v * v;                       // T is now the difference at k + 1
        if (logscale) {
            lpk += lratio + std::log((kd + r) / (kd + 1.0));
            pk = std::exp(lpk);
            // leave log scale only once pk is comfortably NORMAL: the first
            // nonzero exp(lpk) is a subnormal with almost no significand, and
            // seeding the multiplicative recurrence there was measured one
            // family over to carry a 2.5x error to the mode
            if (lpk > -640.0) logscale = false;
        } else {
            pk *= (kd + r) / (kd + 1.0) * ratio;
        }
    }
    return (cum > 0) ? s / cum : 0.0;   // the difference at k = 0 is zero
}

inline void nb1_G_derivs1(double mu, double th, double* G) {
    double r = mu / th;
    double ratio = th / (1.0 + th);
    double lratio = std::log(th) - std::log1p(th);
    double L = std::log1p(th), op = 1.0 + th;
    double cap = 100.0 + mu + 20.0 * d7::sqrt_cr(mu * (1.0 + th))
                 + 40.0 / (-std::log(ratio));
    int kmax = (int) std::min(cap, 2.0e9);
    double lpk = -r * std::log1p(th);
    bool logscale = (lpk <= -640.0);
    double pk = std::exp(lpk);
    double S = 0.0, T = 0.0, U = 0.0, cum = 0.0;
    double s[3] = {0, 0, 0};
    for (int k = 0; k <= kmax; ++k) {
        double kd = (double) k;
        double a = S - L;
        double b = (kd - mu) / (th * op);
        s[0] += pk * T;
        s[1] += pk * (a * T + U);
        s[2] += pk * b * T;
        cum += pk;
        if (cum >= 1.0 - 1e-12) break;
        double v = 1.0 / (r + kd), v2 = v * v;
        S += v;
        T -= v2;
        U += 2.0 * v2 * v;
        if (logscale) {
            lpk += lratio + std::log((kd + r) / (kd + 1.0));
            pk = std::exp(lpk);
            if (lpk > -640.0) logscale = false;
        } else {
            pk *= (kd + r) / (kd + 1.0) * ratio;
        }
    }
    for (int j = 0; j < 3; ++j) G[j] = (cum > 0) ? s[j] / cum : 0.0;
}

namespace d7 {

inline double negbin1_score_mu(const NB1parts& z) {
    return z.P * z.rm;
}

inline double negbin1_score_theta(const NB1parts& z) {
    return z.P * z.rt + z.Q;
}

inline double negbin1_hess_mu_mu(const NB1parts& z) {
    return z.Pr * z.rm * z.rm;
}

inline double negbin1_hess_theta_theta(const NB1parts& z) {
    return z.Pr * z.rt * z.rt + 2.0 * z.Pth * z.rt + z.P * z.rtt + z.Qth;
}

// E[l] at one (mu, theta), EPr = E[psi'(Y + r) - psi'(r)] from nb1_E_Pr()
inline double negbin1_expected_mu_mu(double m, double t, double EPr) {
    double rm = 1.0 / t;
    return EPr * rm * rm;
}

inline double negbin1_expected_theta_theta(double m, double t, double EPr) {
    double r = m / t, om = 1.0 + t;
    double rt = -m / (t * t);
    double Pth = -1.0 / om;
    double EQth = r / (om * om) - m / (t * t) + m / (om * om);
    return EPr * rt * rt + 2.0 * Pth * rt + EQth;
}

// the derivatives of E_mm = G/t^2 in mu and of E_tt = m^2 G/t^4 + m h(t) in
// theta, from g = (G, G_r, G_theta) of nb1_G_derivs1(): d_mu = d_r/t and
// d_theta|mu = d_theta|r - (r/t) d_r
inline double negbin1_dexpected_mu_mu_mu(double m, double t, const double* g) {
    double G = g[0], it = 1.0 / t;
    double Gm = g[1] * it;
    double it2 = it * it;
    double a0[3] = {it2, 0.0, 0.0};
    double z0[3] = {0, 0, 0};
    return z0[1] + a0[1] * G + a0[0] * Gm;
}

inline double negbin1_dexpected_theta_theta_theta(double m, double t,
                                                  const double* g) {
    double r = m / t;
    double G = g[0], it = 1.0 / t, rt = r * it;
    double Gt = g[2] - rt * g[1];
    double it2 = it * it, it3 = it2 * it, it4 = it2 * it2, it5 = it4 * it;
    double op = 1.0 + t, iop = 1.0 / op, iop2 = iop * iop, iop3 = iop2 * iop;
    double a1[3] = {m * m * it4, 2.0 * m * it4, -4.0 * m * m * it5};
    double h1 = 2.0 * (-2.0 * it3 * iop - it2 * iop2)
                - it2 * iop2 - 2.0 * it * iop3 + 2.0 * it3 - 2.0 * iop3;
    return m * h1 + a1[2] * G + a1[0] * Gt;
}

inline void negbin1_score_curv(int k, double y, const double* th,
                               double* out) {
    const NB1parts z = nb1_parts(y, th[0], th[1], 2);
    if (k == 0) {
        out[0] = negbin1_score_mu(z);
        out[1] = negbin1_hess_mu_mu(z);
    } else {
        out[0] = negbin1_score_theta(z);
        out[1] = negbin1_hess_theta_theta(z);
    }
}

inline void negbin1_info_dinfo(int k, double y, const double* th,
                               double* out) {
    double m = th[0], t = th[1];
    double EPr = nb1_E_Pr(m, t);
    double g[3];
    nb1_G_derivs1(m, t, g);
    if (k == 0) {
        out[0] = negbin1_expected_mu_mu(m, t, EPr);
        out[1] = negbin1_dexpected_mu_mu_mu(m, t, g);
    } else {
        out[0] = negbin1_expected_theta_theta(m, t, EPr);
        out[1] = negbin1_dexpected_theta_theta_theta(m, t, g);
    }
}

} // namespace d7

#endif
