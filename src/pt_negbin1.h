#ifndef D7_PT_NEGBIN1_H
#define D7_PT_NEGBIN1_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "psi_diff.h"
#include "pt_sqrt.h"
#include "pt_negbin1_asym.h"

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
// W still cancels as theta -> 0: v_j is of order theta and changes sign, and
// W/theta^2 tends to -1/2, so the relative error of E[l_theta,theta] and
// E[l_mu,theta] grew like eps/theta (1e-7 at theta = 1e-8, 1.8e-4 at 1e-12,
// mu = 4) and d_theta E[l_theta,theta] was -8.59 against 1.125 at 1e-8.
// Below nb1_expected_cut() the sums therefore carry
//   G = E[l_theta,theta] = -E[D_Y] - mu q'(theta),  D_k = sum_{j<k} j^2/u_j^2,
// the expectation of h_theta,theta below (E[Y - mu] = 0), whose summands are
// all positive; the cancellation left is between E[D] (about mu/3 + 1/2) and
// mu q' (about -mu/3), a factor mu and not 1/theta, and the D sums carry the
// rounding of about mu + 40 sd terms. Above the cut G = W/theta^2 from the W
// sums. The cut is where the two meet: against 50-digit sums on a grid of mu
// from 4 to 500 and theta from 1e-4 to 0.9 the W route reads 3e-11 and the D
// route 1.5e-13 at mu = 500, theta = 0.01, and 2.4e-12 against 2e-11 at
// theta = 0.05 (vg_T3_nb1ref5.py). Either way the sums return G with its
// derivatives, and E[l_mu,theta] = -theta G/mu.
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
    // capped at 1e6 terms, as negbin2's sums are (negbin2_series_cap()):
    // past it the mass sits where no estimate reads the information -- a
    // joint fit's trial point reached mu = 1.8e308 (the log link saturated)
    // with theta = 1.9e-77, and that one row cost 85 s at the cap of 2e9
    w.kmax = (int) std::min(cap, 1.0e6);
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
            // q' lags q by one power: stop on the power q' adds next
            if (prev < 1e-18) break;
        }
        *q = s;
        *qp = sp;
    } else {
        const double tml = t - std::log1p(t);
        *q = tml / (t * t);
        *qp = 1.0 / (t * (1.0 + t)) - 2.0 * tml / (t * t * t);
    }
}

// q'' and q''' for the derivatives of G below, by
//   q''  = sum_{n>=4} (-1)^n theta^(n-4) (n - 2)(n - 3)/n,
//   q''' = sum_{n>=5} (-1)^n theta^(n-5) (n - 2)(n - 3)(n - 4)/n
// below theta = 0.5, and above it by the closed forms, with T = theta -
// log1p(theta) and s = theta (1 + theta),
//   q''  = -(1 + 2 theta)/s^2 - 2/(theta^2 (1 + theta)) + 6 T/theta^4,
//   q''' = 2 (1 + 2 theta)^2/s^3 + 10/(theta^3 (1 + theta)) - 24 T/theta^5,
// which lose a factor of a few hundred at theta = 0.5 and less above.
inline void nb1_q23(double t, double* qpp, double* qppp) {
    if (t < 0.5) {
        double s2 = 0.0, s3 = 0.0, pw = 1.0, prev = 0.0;   // pw = t^(n-4)
        for (int n = 4; n < 200; ++n) {
            const double sg = (n % 2 == 0) ? 1.0 : -1.0;
            s2 += sg * pw * (n - 2) * (n - 3) / n;
            if (n >= 5) s3 += sg * prev * (n - 2) * (n - 3) * (n - 4) / n;
            prev = pw;
            pw *= t;
            // q''' lags q'' by one power: stop on the power q''' adds next
            if (prev * n * n * n < 1e-18) break;
        }
        *qpp = s2;
        *qppp = s3;
    } else {
        const double T = t - std::log1p(t), s = t * (1.0 + t), w = 1.0 + 2.0 * t;
        const double t2 = t * t, t3 = t2 * t;
        *qpp = -w / (s * s) - 2.0 / (t2 * (1.0 + t)) + 6.0 * T / (t2 * t2);
        *qppp = 2.0 * w * w / (s * s * s) + 10.0 / (t3 * (1.0 + t))
                - 24.0 * T / (t3 * t2);
    }
}

// The dispersion below which the scores and second derivatives are read
// from the sums in mu + i theta (nb1_low_exact) rather than from
// nb1_parts(): the cut of the higher orders, nb1_exact_cut() in
// R/negbin1_higher.R. Above it the two agree to 1e-11 or better.
inline double nb1_low_cut() { return 1.0; }

// The dispersion below which the expected information and its derivatives
// are summed in G = -E[D] - mu q' (nb1_dsums0/1/2) rather than in W.
inline double nb1_expected_cut() { return 0.02; }

// The scores and second derivatives from the form in mu + i theta, by the
// formulas for s_mu, s_theta, h_mu,mu, h_mu,theta and h_theta,theta above
// nb1_q(). No term of order mu/theta^2 is formed, so every component keeps
// its accuracy down to the Poisson limit; the cost is y terms an
// observation. out holds s_mu, s_theta and, at order 2, h_mu,mu,
// h_mu,theta, h_theta,theta.
inline void nb1_low_exact(double y, double mu, double th, int order,
                          double* out) {
    double B = 0.0, a = 0.0, C = 0.0, D = 0.0;
    const long k = (long) y;
    if (k >= 1) a = 1.0 / (mu * mu);
    for (long i = 1; i < k; ++i) {
        const double di = (double) i;
        const double iu = 1.0 / (mu + di * th);
        B += di * iu;
        if (order >= 2) {
            const double iu2 = iu * iu;
            a += iu2;
            C += di * iu2;
            D += di * di * iu2;
        }
    }
    double q, qp;
    nb1_q(th, &q, &qp);
    const double om = 1.0 + th, km = y - mu;
    out[0] = km / mu - (th / mu) * B + th * q;
    out[1] = -km / om + B - mu * q;
    if (order < 2) return;
    out[2] = -a;
    out[3] = 1.0 / om - C - q;
    out[4] = km / (om * om) - D - mu * qp;
}

// A and W at one (mu, theta)
inline void nb1_wsums0(double mu, double th, double* A, double* W) {
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
inline void nb1_wsums1(double mu, double th, double* A, double* W) {
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
inline void nb1_wsums2(double mu, double th, double* A, double* W) {
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

// A and G = E[l_theta,theta] below nb1_expected_cut(), from the D sums
inline void nb1_dsums0(double mu, double th, double* A, double* G) {
    NB1mass w = nb1_mass(mu, th);
    const double r = w.r, ratio = w.ratio, lratio = w.lratio, mode = w.mode;
    double pk = w.pk, lpk = w.lpk;
    bool logscale = w.logscale;
    double q, qp;
    nb1_q(th, &q, &qp);
    double a = 0.0, Dv = 0.0, sa = 0.0, sd = 0.0;
    for (int k = 0; k <= w.kmax; ++k) {
        const double kd = (double) k;
        sa += pk * a;
        sd += pk * Dv;
        if (negbin1_series_done(kd, pk, mode)) break;
        const double u = mu + th * kd, iu = 1.0 / u, iu2 = iu * iu;
        a += iu2;
        Dv += kd * kd * iu2;
        nb1_mass_step(kd, r, ratio, lratio, pk, lpk, logscale);
    }
    *A = sa;
    *G = -sd - mu * qp;
}

// the same with the first derivatives, each as (value, d_mu, d_theta)
inline void nb1_dsums1(double mu, double th, double* A, double* G) {
    NB1mass w = nb1_mass(mu, th);
    const double r = w.r, ratio = w.ratio, lratio = w.lratio, mode = w.mode;
    double pk = w.pk, lpk = w.lpk;
    bool logscale = w.logscale;
    const double iop = 1.0 / (1.0 + th), imu = 1.0 / mu, thmu = th * imu;
    double q, qp, qpp, qppp;
    nb1_q(th, &q, &qp);
    nb1_q23(th, &qpp, &qppp);
    const double thq = th * q, muq = mu * q;
    double a = 0.0, am = 0.0, at = 0.0, Dv = 0.0, Dm = 0.0, Dt = 0.0, B = 0.0;
    double sa[3] = {0, 0, 0}, sd[3] = {0, 0, 0};
    for (int k = 0; k <= w.kmax; ++k) {
        const double kd = (double) k, p = pk;
        const double sm = (kd - mu) * imu - thmu * B + thq;
        const double st = -(kd - mu) * iop + B - muq;
        sa[0] += p * a;
        sa[1] += p * (am + a * sm);
        sa[2] += p * (at + a * st);
        sd[0] += p * Dv;
        sd[1] += p * (Dm + Dv * sm);
        sd[2] += p * (Dt + Dv * st);
        if (negbin1_series_done(kd, p, mode)) break;
        const double u = mu + th * kd, iu = 1.0 / u, iu2 = iu * iu, iu3 = iu2 * iu;
        const double k2 = kd * kd;
        a += iu2;
        am -= 2.0 * iu3;
        at -= 2.0 * kd * iu3;
        Dv += k2 * iu2;
        Dm -= 2.0 * k2 * iu3;
        Dt -= 2.0 * k2 * kd * iu3;
        B += kd * iu;
        nb1_mass_step(kd, r, ratio, lratio, pk, lpk, logscale);
    }
    for (int j = 0; j < 3; ++j) A[j] = sa[j];
    G[0] = -sd[0] - mu * qp;
    G[1] = -sd[1] - qp;
    G[2] = -sd[2] - mu * qpp;
}

// the same with the first and second derivatives, each as (value, d_mu,
// d_theta, d_mu,mu, d_mu,theta, d_theta,theta)
inline void nb1_dsums2(double mu, double th, double* A, double* G) {
    NB1mass w = nb1_mass(mu, th);
    const double r = w.r, ratio = w.ratio, lratio = w.lratio, mode = w.mode;
    double pk = w.pk, lpk = w.lpk;
    bool logscale = w.logscale;
    const double iop = 1.0 / (1.0 + th), iop2 = iop * iop;
    const double imu = 1.0 / mu, thmu = th * imu;
    double q, qp, qpp, qppp;
    nb1_q(th, &q, &qp);
    nb1_q23(th, &qpp, &qppp);
    const double thq = th * q, muq = mu * q, muqp = mu * qp;
    double f[2][6] = {{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};
    double s[2][6] = {{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};
    double B = 0.0, C = 0.0;
    for (int k = 0; k <= w.kmax; ++k) {
        const double kd = (double) k, p = pk;
        const double sm = (kd - mu) * imu - thmu * B + thq;
        const double st = -(kd - mu) * iop + B - muq;
        const double hmm = -f[0][0];
        const double hmt = iop - C - q;
        const double htt = (kd - mu) * iop2 - f[1][0] - muqp;
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
        const double k2 = kd * kd;
        double* a = f[0];
        a[0] += iu2;
        a[1] -= 2.0 * iu3;
        a[2] -= 2.0 * kd * iu3;
        a[3] += 6.0 * iu4;
        a[4] += 6.0 * kd * iu4;
        a[5] += 6.0 * k2 * iu4;
        double* D = f[1];
        D[0] += k2 * iu2;
        D[1] -= 2.0 * k2 * iu3;
        D[2] -= 2.0 * k2 * kd * iu3;
        D[3] += 6.0 * k2 * iu4;
        D[4] += 6.0 * k2 * kd * iu4;
        D[5] += 6.0 * k2 * k2 * iu4;
        B += kd * iu;
        C += kd * iu2;
        nb1_mass_step(kd, r, ratio, lratio, pk, lpk, logscale);
    }
    for (int j = 0; j < 6; ++j) A[j] = s[0][j];
    G[0] = -s[1][0] - mu * qp;
    G[1] = -s[1][1] - qp;
    G[2] = -s[1][2] - mu * qpp;
    G[3] = -s[1][3];
    G[4] = -s[1][4] - qpp;
    G[5] = -s[1][5] - mu * qppp;
}

// A and G = E[l_theta,theta] at one (mu, theta), from the D sums below
// nb1_expected_cut() and from G = W/theta^2 above it
// LARGE THETA. At fixed mu the mass has about theta terms and a ratio
// theta/(1 + theta) close to one. The series of pt_negbin1_asym.h replaces the
// sums wherever it meets 1e-13 (nb1_asym_ok(): theta >= 2, mu/theta <= 0.3),
// where it costs about 0.5 us against 18-30 us for the sums at theta = 50
// and none of their length at large theta (11 s a row at theta = 1e8).
inline bool nb1_use_asym(double mu, double th) {
    return nb1_asym_ok(mu, th);
}

inline void nb1_sums0(double mu, double th, double* A, double* G) {
    if (nb1_use_asym(mu, th)) { nb1_asym_sums0(mu, th, A, G); return; }
    if (th < nb1_expected_cut()) { nb1_dsums0(mu, th, A, G); return; }
    double W;
    nb1_wsums0(mu, th, A, &W);
    *G = W / (th * th);
}

inline void nb1_sums1(double mu, double th, double* A, double* G) {
    if (nb1_use_asym(mu, th)) { nb1_asym_sums1(mu, th, A, G); return; }
    if (th < nb1_expected_cut()) { nb1_dsums1(mu, th, A, G); return; }
    double W[3];
    nb1_wsums1(mu, th, A, W);
    const double i2 = 1.0 / (th * th), i3 = i2 / th;
    G[0] = W[0] * i2;
    G[1] = W[1] * i2;
    G[2] = W[2] * i2 - 2.0 * W[0] * i3;
}

inline void nb1_sums2(double mu, double th, double* A, double* G) {
    if (nb1_use_asym(mu, th)) { nb1_asym_sums2(mu, th, A, G); return; }
    if (th < nb1_expected_cut()) { nb1_dsums2(mu, th, A, G); return; }
    double W[6];
    nb1_wsums2(mu, th, A, W);
    const double i2 = 1.0 / (th * th), i3 = i2 / th, i4 = i2 * i2;
    G[0] = W[0] * i2;
    G[1] = W[1] * i2;
    G[2] = W[2] * i2 - 2.0 * W[0] * i3;
    G[3] = W[3] * i2;
    G[4] = W[4] * i2 - 2.0 * W[1] * i3;
    G[5] = W[5] * i2 - 4.0 * W[2] * i3 + 6.0 * W[0] * i4;
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
// E[l_theta,theta] = G, E[l_mu,theta] = -theta G/mu
inline double negbin1_expected_mu_mu(double A) {
    return -A;
}

inline double negbin1_expected_theta_theta(double t, double G) {
    return G;
}

inline double negbin1_expected_mu_theta(double m, double t, double G) {
    return -t * G / m;
}

// Its first derivatives, from nb1_sums1()'s (value, d_mu, d_theta) of A
// and G
inline double negbin1_dexpected_mu_mu_mu(const double* A) {
    return -A[1];
}

inline double negbin1_dexpected_mu_mu_theta(const double* A) {
    return -A[2];
}

inline double negbin1_dexpected_theta_theta_mu(double t, const double* G) {
    return G[1];
}

inline double negbin1_dexpected_theta_theta_theta(double t, const double* G) {
    return G[2];
}

inline double negbin1_dexpected_mu_theta_mu(double m, double t,
                                            const double* G) {
    return -t * G[1] / m + t * G[0] / (m * m);
}

inline double negbin1_dexpected_mu_theta_theta(double m, double t,
                                               const double* G) {
    return -G[0] / m - t * G[2] / m;
}

// The scores and second derivatives at one observation, by the route that
// is accurate at the dispersion given: below nb1_low_cut() the sums in
// mu + i theta, above it nb1_parts(). A response that is not a count goes
// to nb1_parts(), as before.
inline void negbin1_low_orders(double y, double m, double t, int order,
                               double* out) {
    if (t < nb1_low_cut() && y >= 0.0 && y == std::floor(y)) {
        nb1_low_exact(y, m, t, order, out);
        return;
    }
    const NB1parts z = nb1_parts(y, m, t, order);
    out[0] = negbin1_score_mu(z);
    out[1] = negbin1_score_theta(z);
    if (order < 2) return;
    out[2] = negbin1_hess_mu_mu(z);
    out[3] = (z.Pr * z.rt + z.Pth) * z.rm + z.P * z.rmt;
    out[4] = negbin1_hess_theta_theta(z);
}

inline void negbin1_score_curv(int k, double y, const double* th,
                               double* out) {
    double v[5];
    negbin1_low_orders(y, th[0], th[1], 2, v);
    if (k == 0) {
        out[0] = v[0];
        out[1] = v[2];
    } else {
        out[0] = v[1];
        out[1] = v[4];
    }
}

inline void negbin1_info_dinfo(int k, double y, const double* th,
                               double* out) {
    double m = th[0], t = th[1];
    double A[3], G[3];
    nb1_sums1(m, t, A, G);
    if (k == 0) {
        out[0] = negbin1_expected_mu_mu(A[0]);
        out[1] = negbin1_dexpected_mu_mu_mu(A);
    } else {
        out[0] = negbin1_expected_theta_theta(t, G[0]);
        out[1] = negbin1_dexpected_theta_theta_theta(t, G);
    }
}

// The routers of the joint registry, one per order, for (mu, theta) in that
// order. Each fills a full P x P matrix (P = 2) column-major, out[a + P b],
// and the derivative of the expected information a P x P x P array, out[a +
// P b + P P c] being the derivative of E[l_ab] in parameter c. The grad and
// the hess call negbin1_low_orders(), as negbin1.cpp does, so the route to
// the exact sums (integer y below theta = 1) and to nb1_parts() is the same.
inline void negbin1_joint_grad(double y, const double* th, double* out) {
    double v[5];
    negbin1_low_orders(y, th[0], th[1], 1, v);
    out[0] = v[0];
    out[1] = v[1];
}

inline void negbin1_joint_hess(double y, const double* th, double* out) {
    double v[5];
    negbin1_low_orders(y, th[0], th[1], 2, v);
    out[0] = v[2];
    out[3] = v[4];
    out[1] = out[2] = v[3];
}

// The expected information from nb1_sums0(), as nb1_E_parts() in
// negbin1.cpp forms it.
inline void negbin1_joint_info(double y, const double* th, double* out) {
    double m = th[0], t = th[1];
    double A, G;
    nb1_sums0(m, t, &A, &G);
    out[0] = negbin1_expected_mu_mu(A);
    out[3] = negbin1_expected_theta_theta(t, G);
    out[1] = out[2] = negbin1_expected_mu_theta(m, t, G);
}

// The derivatives of the expected information from nb1_sums1(), in the
// order of negbin1_dexpected1_cpp().
inline void negbin1_joint_dinfo(double y, const double* th, double* out) {
    double m = th[0], t = th[1];
    double A[3], G[3];
    nb1_sums1(m, t, A, G);
    for (int j = 0; j < 8; ++j) out[j] = 0.0;
    out[0] = negbin1_dexpected_mu_mu_mu(A);
    out[4] = negbin1_dexpected_mu_mu_theta(A);
    out[2] = out[1] = negbin1_dexpected_mu_theta_mu(m, t, G);
    out[6] = out[5] = negbin1_dexpected_mu_theta_theta(m, t, G);
    out[3] = negbin1_dexpected_theta_theta_mu(t, G);
    out[7] = negbin1_dexpected_theta_theta_theta(t, G);
}

} // namespace d7

#endif
