#ifndef D7_PT_NEGBIN1_H
#define D7_PT_NEGBIN1_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "psi_diff.h"
#include "pt_sqrt.h"

// The negative binomial with a variance linear in the mean, in (mu, theta),
// r = mu/theta: the parts of the log-mass (nb1_parts), the support sums that
// carry its expected information (nb1_sums0, nb1_sums1, nb1_sums2), and one
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

// THE EXPECTED INFORMATION AND ITS DERIVATIVES, from two sums over the
// support written in (mu, theta) themselves. With u_j = mu + theta j,
//
//   A = E[a_Y],  a_k = sum_{j<k} 1/u_j^2,
//   W = E[V_Y],  V_k = sum_{j<k} v_j,
//   v_j = 1/(1 + theta) - mu^2/u_j^2
//       = theta (2 mu j + theta j^2 - mu^2) / ((1 + theta) u_j^2),
//
// the three components are E[l_mu,mu] = -A, E[l_theta,theta] = W/theta^2
// and E[l_mu,theta] = -W/(mu theta). The first is E[psi'(Y + r) -
// psi'(r)]/theta^2 by the recurrence psi'(r + k) - psi'(r) = -sum_{j<k}
// 1/(r + j)^2, r = mu/theta. The other two were E[psi'] composed with closed
// forms, m^2 G/theta^4 + m/(theta^2 (1 + theta)), whose pieces are of order
// mu/theta^2 and whose sum is of order one; the closed form is itself an
// expectation, mu/(theta^2 (1 + theta)) = E[Y]/(theta^2 (1 + theta)), and
// the two merge into the one summand v_j, an exact polynomial over u_j^2.
//
// The derivatives of A and W are sums over the same mass with the score of
// the mass beside the summand,
//   d_c F   = E[f_c + f s_c],
//   d_ce F  = E[f_ce + f_c s_e + f_e s_c + f (s_c s_e + h_ce)],
// taken in (mu, theta) directly: the earlier kernels took them in
// (r, theta) and changed coordinates by d_theta|mu = d_theta|r - (r/theta)
// d_r, whose two terms are of order mu/theta^2 apiece. Against sums at 50
// digits on 16 cells (mu 0.5 to 100, theta 0.05 to 20) the worst error
// relative to the largest component of its order went from 7.4e-7
// (d_theta^2 E[l_theta,theta] at mu = 100, theta = 0.05) to 4.3e-10, and
// d_mu E[l_theta,theta] there, a component of size 2e-6, from 4.5e-6 of
// itself to 2e-10.
//
// The scores and second derivatives of the log-mass, written so that no
// term of order mu/theta^2 is formed, with B_k = sum_{j<k} j/u_j,
// C_k = sum_{j<k} j/u_j^2, D_k = sum_{j<k} j^2/u_j^2 and
// q = (theta - log1p(theta))/theta^2:
//   s_mu        = (k - mu)/mu - (theta/mu) B_k + theta q,
//   s_theta     = -(k - mu)/(1 + theta) + B_k - mu q,
//   h_mu,mu     = -a_k,
//   h_mu,theta  = 1/(1 + theta) - C_k - q,
//   h_theta,th  = (k - mu)/(1 + theta)^2 - D_k - mu q'.
// q and q' are their power series below theta = 0.25 (nb1_q()), where the
// closed forms lose a factor 1/theta and 1/theta^2.
//
// The mass comes from the recurrence p_{k+1} = p_k (k + r)/(k + 1) *
// theta/(1 + theta), carried in log scale until it is representable, and the
// loop stops by negbin1_series_done(). An earlier version located the tail
// by a far quantile, and the call was R::qnbinom, which is off limits here for the reason
// nb_E_trigamma states one family over: these helpers run inside d7::par_for
// workers, and qnbinom's search reaches pbeta, whose warning path calls into
// the R API and killed the process from a worker thread on four of five CI
// platforms.
//
// The hard cap covers the geometric tail, whose decay ratio is
// theta/(1 + theta); the loop almost always breaks on the mass long before.
// THE STOPPING RULE OF THE SERIES BELOW, negbin2's (pt_negbin2.h): past the
// mode, and not before k = 100, once p_k (1 + k)^2 falls below 1e-17. The
// rule it replaces stopped at 1 - 1e-12 of the accumulated mass; against
// 50-digit sums on 16 cells (mu 0.5 to 100, theta 0.05 to 20) that left up
// to 1.5e-2 on d_mu E[l_theta,theta] at mu = 100, theta = 0.05, and 7.8e-5
// at mu = 10. The mode of the mass is (r - 1) theta for r = mu/theta above
// one.
inline bool negbin1_series_done(double kd, double pk, double mode) {
    return kd >= 100.0 && kd > mode && pk * (1.0 + kd) * (1.0 + kd) <= 1e-17;
}

inline double negbin1_series_mode(double r, double th) {
    return (r > 1.0) ? (r - 1.0) * th : 0.0;
}

// the mass at k = 0 and the recurrence's constants, shared by the three sums
struct NB1mass {
    double r, ratio, lratio, mode, lpk, pk;
    int kmax;
    bool logscale;
};

inline NB1mass nb1_mass(double mu, double th) {
    NB1mass w;
    w.r = mu / th;
    w.ratio = th / (1.0 + th);
    w.lratio = std::log(th) - std::log1p(th);
    double cap = 100.0 + mu + 20.0 * d7::sqrt_cr(mu * (1.0 + th))
                 + 80.0 / (-std::log(w.ratio));
    w.kmax = (int) std::min(cap, 2.0e9);
    w.mode = negbin1_series_mode(w.r, th);
    // The log scale is left on the SAME threshold the loop switches at, and
    // not on exp(lpk) merely being nonzero: at mu = 918.832, theta = 0.5 the
    // size is 1837.7 and lpk is -745.1, so P(Y = 0) is a subnormal with
    // almost no significand. Seeding the multiplicative recurrence there put
    // mu_theta at -2.73e-02 where it is 1.2089e-04.
    w.lpk = -w.r * std::log1p(th);
    w.logscale = (w.lpk <= -640.0);
    w.pk = std::exp(w.lpk);
    return w;
}

// the step from p_k to p_{k+1}, on the caller's locals
inline void nb1_mass_step(double kd, double r, double ratio, double lratio,
                          double& pk, double& lpk, bool& logscale) {
    if (logscale) {
        lpk += lratio + std::log((kd + r) / (kd + 1.0));
        pk = std::exp(lpk);
        // leave log scale only once pk is comfortably NORMAL: the first
        // nonzero exp(lpk) is a subnormal with almost no significand
        if (lpk > -640.0) logscale = false;
    } else {
        pk *= (kd + r) / (kd + 1.0) * ratio;
    }
}

// q = (theta - log1p(theta))/theta^2 and its derivative q' = 1/(theta (1 +
// theta)) - 2 (theta - log1p(theta))/theta^3, by
//   q  = sum_{n>=2} (-1)^n theta^(n-2)/n,
//   q' = sum_{n>=3} (-1)^n theta^(n-3) (n - 2)/n
// below theta = 0.25, where the closed forms lose a factor 1/theta and
// 1/theta^2
inline void nb1_q(double t, double* q, double* qp) {
    if (t < 0.25) {
        double s = 0.0, sp = 0.0, pw = 1.0, prev = 0.0;   // pw = t^(n-2)
        for (int n = 2; n < 80; ++n) {
            const double sg = (n % 2 == 0) ? 1.0 : -1.0;
            s += sg * pw / n;
            if (n >= 3) sp += sg * prev * (n - 2) / n;
            prev = pw;
            pw *= t;
            if (pw < 1e-18) break;
        }
        *q = s;
        *qp = sp;
    } else {
        const double tml = t - std::log1p(t);
        *q = tml / (t * t);
        *qp = 1.0 / (t * (1.0 + t)) - 2.0 * tml / (t * t * t);
    }
}

// A and W at one (mu, theta)
inline void nb1_sums0(double mu, double th, double* A, double* W) {
    NB1mass w = nb1_mass(mu, th);
    const double r = w.r, ratio = w.ratio, lratio = w.lratio, mode = w.mode;
    double pk = w.pk, lpk = w.lpk;
    bool logscale = w.logscale;
    const double thop = th * (1.0 / (1.0 + th)), mu2 = mu * mu;
    double a = 0.0, V = 0.0, sa = 0.0, sv = 0.0;
    for (int k = 0; k <= w.kmax; ++k) {
        const double kd = (double) k;
        sa += pk * a;
        sv += pk * V;
        if (negbin1_series_done(kd, pk, mode)) break;
        const double u = mu + th * kd, iu = 1.0 / u, iu2 = iu * iu;
        a += iu2;
        V += thop * (2.0 * mu * kd + th * kd * kd - mu2) * iu2;
        nb1_mass_step(kd, r, ratio, lratio, pk, lpk, logscale);
    }
    *A = sa;
    *W = sv;
}

// A and W with their first derivatives, each as (value, d_mu, d_theta)
inline void nb1_sums1(double mu, double th, double* A, double* W) {
    NB1mass w = nb1_mass(mu, th);
    const double r = w.r, ratio = w.ratio, lratio = w.lratio, mode = w.mode;
    double pk = w.pk, lpk = w.lpk;
    bool logscale = w.logscale;
    const double op = 1.0 + th, iop = 1.0 / op, iop2 = iop * iop;
    const double thop = th * iop, mu2 = mu * mu, imu = 1.0 / mu, thmu = th * imu;
    double q, qp;
    nb1_q(th, &q, &qp);
    const double thq = th * q, muq = mu * q;
    double a = 0.0, am = 0.0, at = 0.0, V = 0.0, Vm = 0.0, Vt = 0.0, B = 0.0;
    double sa[3] = {0, 0, 0}, sv[3] = {0, 0, 0};
    for (int k = 0; k <= w.kmax; ++k) {
        const double kd = (double) k, p = pk;
        const double sm = (kd - mu) * imu - thmu * B + thq;
        const double st = -(kd - mu) * iop + B - muq;
        sa[0] += p * a;
        sa[1] += p * (am + a * sm);
        sa[2] += p * (at + a * st);
        sv[0] += p * V;
        sv[1] += p * (Vm + V * sm);
        sv[2] += p * (Vt + V * st);
        if (negbin1_series_done(kd, p, mode)) break;
        const double u = mu + th * kd, iu = 1.0 / u, iu2 = iu * iu, iu3 = iu2 * iu;
        a += iu2;
        am -= 2.0 * iu3;
        at -= 2.0 * kd * iu3;
        V += thop * (2.0 * mu * kd + th * kd * kd - mu2) * iu2;
        Vm -= 2.0 * mu * th * kd * iu3;
        Vt += 2.0 * mu2 * kd * iu3 - iop2;
        B += kd * iu;
        nb1_mass_step(kd, r, ratio, lratio, pk, lpk, logscale);
    }
    for (int j = 0; j < 3; ++j) { A[j] = sa[j]; W[j] = sv[j]; }
}

// A and W with their first and second derivatives, each as (value, d_mu,
// d_theta, d_mu,mu, d_mu,theta, d_theta,theta)
inline void nb1_sums2(double mu, double th, double* A, double* W) {
    NB1mass w = nb1_mass(mu, th);
    const double r = w.r, ratio = w.ratio, lratio = w.lratio, mode = w.mode;
    double pk = w.pk, lpk = w.lpk;
    bool logscale = w.logscale;
    const double op = 1.0 + th, iop = 1.0 / op, iop2 = iop * iop,
                 iop3 = iop2 * iop;
    const double thop = th * iop, mu2 = mu * mu, imu = 1.0 / mu, thmu = th * imu;
    double q, qp;
    nb1_q(th, &q, &qp);
    const double thq = th * q, muq = mu * q, muqp = mu * qp;
    double f[2][6] = {{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};
    double s[2][6] = {{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};
    double B = 0.0, C = 0.0, D = 0.0;
    for (int k = 0; k <= w.kmax; ++k) {
        const double kd = (double) k, p = pk;
        const double sm = (kd - mu) * imu - thmu * B + thq;
        const double st = -(kd - mu) * iop + B - muq;
        const double hmm = -f[0][0];
        const double hmt = iop - C - q;
        const double htt = (kd - mu) * iop2 - D - muqp;
        for (int e = 0; e < 2; ++e) {
            const double* g = f[e];
            double* o = s[e];
            o[0] += p * g[0];
            o[1] += p * (g[1] + g[0] * sm);
            o[2] += p * (g[2] + g[0] * st);
            o[3] += p * (g[3] + 2.0 * g[1] * sm + g[0] * (sm * sm + hmm));
            o[4] += p * (g[4] + g[1] * st + g[2] * sm + g[0] * (sm * st + hmt));
            o[5] += p * (g[5] + 2.0 * g[2] * st + g[0] * (st * st + htt));
        }
        if (negbin1_series_done(kd, p, mode)) break;
        const double u = mu + th * kd, iu = 1.0 / u, iu2 = iu * iu,
                     iu3 = iu2 * iu, iu4 = iu2 * iu2;
        double* a = f[0];
        a[0] += iu2;
        a[1] -= 2.0 * iu3;
        a[2] -= 2.0 * kd * iu3;
        a[3] += 6.0 * iu4;
        a[4] += 6.0 * kd * iu4;
        a[5] += 6.0 * kd * kd * iu4;
        double* V = f[1];
        V[0] += thop * (2.0 * mu * kd + th * kd * kd - mu2) * iu2;
        V[1] -= 2.0 * mu * th * kd * iu3;
        V[2] += 2.0 * mu2 * kd * iu3 - iop2;
        V[3] -= 2.0 * th * kd * (th * kd - 2.0 * mu) * iu4;
        V[4] -= 2.0 * mu * kd * (mu - 2.0 * th * kd) * iu4;
        V[5] += 2.0 * iop3 - 6.0 * mu2 * kd * kd * iu4;
        B += kd * iu;
        C += kd * iu2;
        D += kd * kd * iu2;
        nb1_mass_step(kd, r, ratio, lratio, pk, lpk, logscale);
    }
    for (int j = 0; j < 6; ++j) { A[j] = s[0][j]; W[j] = s[1][j]; }
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

// The expected Hessian from nb1_sums0(): E[l_mu,mu] = -A,
// E[l_theta,theta] = W/theta^2, E[l_mu,theta] = -W/(mu theta)
inline double negbin1_expected_mu_mu(double A) {
    return -A;
}

inline double negbin1_expected_theta_theta(double t, double W) {
    return W / (t * t);
}

inline double negbin1_expected_mu_theta(double m, double t, double W) {
    return -W / (m * t);
}

// Its first derivatives, from nb1_sums1()'s (value, d_mu, d_theta)
inline double negbin1_dexpected_mu_mu_mu(const double* A) {
    return -A[1];
}

inline double negbin1_dexpected_mu_mu_theta(const double* A) {
    return -A[2];
}

inline double negbin1_dexpected_theta_theta_mu(double t, const double* W) {
    return W[1] / (t * t);
}

inline double negbin1_dexpected_theta_theta_theta(double t, const double* W) {
    return W[2] / (t * t) - 2.0 * W[0] / (t * t * t);
}

inline double negbin1_dexpected_mu_theta_mu(double m, double t,
                                            const double* W) {
    return -W[1] / (m * t) + W[0] / (m * m * t);
}

inline double negbin1_dexpected_mu_theta_theta(double m, double t,
                                               const double* W) {
    return -W[2] / (m * t) + W[0] / (m * t * t);
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
    double A[3], W[3];
    nb1_sums1(m, t, A, W);
    if (k == 0) {
        out[0] = negbin1_expected_mu_mu(A[0]);
        out[1] = negbin1_dexpected_mu_mu_mu(A);
    } else {
        out[0] = negbin1_expected_theta_theta(t, W[0]);
        out[1] = negbin1_dexpected_theta_theta_theta(t, W);
    }
}

} // namespace d7

#endif
