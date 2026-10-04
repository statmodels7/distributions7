#ifndef D7_PT_BETABINOM_H
#define D7_PT_BETABINOM_H

#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <vector>
#include "psi_diff.h"

// The beta-binomial in its two parametrizations: betabinom1 in the mean
// proportion and the dispersion (mu, sigma), betabinom2 in the shapes
// (alpha, beta), both with size N. What the kernels of betabinom.cpp and
// dexpected_kernels.cpp and the scalar registry share: the log-mass, the
// derivatives in the shapes (bb_shape_derivs), the map to (mu, sigma)
// (bb_map), the power sums over the support (BBSums), one function per
// component built on them, and the expected information's diagonal with its
// derivative, summed exactly over the support.

namespace d7 {

// The log-mass, by whichever of two routes is accurate at the shapes given.
// The two beta functions are of magnitude S log S and their difference is of
// order one, so the ordinary route carries an absolute error of eps times
// that magnitude; past 1e-8 the shifts, being integers, are expanded as exact
// sums of logarithms, which form nothing larger than n log S. Without this
// the difference of the two beta functions collapses to zero at large shapes
// and every mass comes back as one.
// R's lchoose() (src/nmath/choose.c, R 4.6.0) at a non-negative integer n
// and an integer k, its branches and expressions without the R_CheckStack()
// it calls on every entry: that check reads the calling thread's stack
// bounds and, from a worker thread, aborts with "C stack usage is too close
// to the limit", whatever the arguments. The values are lchoose()'s.
inline double lchoose_int(double n, double k) {
    if (k < 2) {
        if (k < 0) return R_NegInf;
        if (k == 0) return 0.;
        return std::log(std::fabs(n));
    }
    if (n < k) return R_NegInf;
    if (n - k < 2) return lchoose_int(n, n - k);
    return -std::log(n + 1.) - R::lbeta(n - k + 1., k + 1.);
}

// Every caller reaches this with an integer y in 0..n (the support guard of
// betabinom_logpmf_cpp and of the registry's log-density, or the sums over
// the support), which lchoose_int() requires.
inline double bb_log_mass(double y, double A, double B, double n) {
    double S = A + B;
    static const double eps = std::numeric_limits<double>::epsilon();
    if (R_FINITE(S) && R::lgammafn(S + n) * eps < 1e-8) {
        return lchoose_int(n, y) + (R::lbeta(y + A, n - y + B) - R::lbeta(A, B));
    }
    double s1 = 0.0, s2 = 0.0, s3 = 0.0;
    int N = (int) n;
    for (int j = 0; j < N; j++) {
        if (j < y)          s1 += std::log(A + j);
        if (j < n - y)      s2 += std::log(B + j);
        s3 += std::log(S + j);
    }
    return lchoose_int(n, y) + s1 + s2 - s3;
}

struct BBderiv {
    double lA, lB;              // first order in the shapes
    double lAA, lAB, lBB;       // second order
};

// EACH CANCELLATION WRITTEN OUT.  As the concentration S = A + B grows the
// family tends to the binomial and every derivative in the shapes vanishes,
// so each was a difference of digammas at arguments a whole SIZE apart:
// measured, dl/dA is wrong by 5.6e-05 at S = 1e6, EXACTLY ZERO at 1e9 where
// the value is 1.9e-17, and 1.9e+08 out at 1e12 -- while a fit at a true
// concentration of 3000 reports 1.7e+08 and one on binomial data 3.1e+09.
//
// The log-mass above has carried the exact form since 0.20.0, summing
// log(A+j) over the support rather than differencing two lbeta; these are
// the derivatives of that sum and had not followed it.  With
// psi(x+k) - psi(x) = psi_A_rest(k,x) + log1p(k/x), the two logarithms
// combine into one log1p of a small quantity:
//   log1p(y/A) - log1p(n/S) = log1p((y S - A n)/(A(S+n)))
// and the same for B.  Nothing here is O(n): the sums the log-mass needs are
// replaced by the series, which is O(1) in the size.
// The first order alone, and the second alone: a caller reading one order
// computes no polygamma of the other.
struct BBd1 { double lA, lB; };
struct BBd2 { double lAA, lAB, lBB; };

inline BBd1 bb_shape_d1(double y, double n, double A, double B) {
    double S = A + B;
    double Sn = S + n;
    BBd1 d;
    double aS  = psi_A_rest(n, S);
    d.lA  = psi_A_rest(y, A) - aS +
        std::log1p((y * S - A * n) / (A * Sn));
    d.lB  = psi_A_rest(n - y, B) - aS +
        std::log1p(((n - y) * S - B * n) / (B * Sn));
    return d;
}

// psi(x + k) - psi(x): for an integer shift up to 32 the sum of
// 1/(x + i), i < k, of positive terms; otherwise two digammas below the
// series' crossover of psi_A_rest(), the series plus log1p(k/x) above it
inline double psi_diff1(double k, double x) {
    if (x >= 100.0) return psi_A_rest(k, x) + std::log1p(k / x);
    if (k <= 32.0 && k >= 0.0 && k == std::floor(k)) {
        double s = 0.0;
        for (int i = 0; i < (int) k; ++i) s += 1.0 / (x + i);
        return s;
    }
    return R::digamma(x + k) - R::digamma(x);
}

// bb_shape_d1() with the logarithms combined only where they cancel: a
// shape at or above 100 makes its component a difference of two quantities
// of size n/S, which the combined log1p of bb_shape_d1() keeps; below it the
// component is a difference of digammas, which does not cancel, and the
// log1p pair adds nothing but cost
inline BBd1 bb_shape_d1_mixed(double y, double n, double A, double B) {
    if (A >= 100.0 && B >= 100.0) return bb_shape_d1(y, n, A, B);
    double S = A + B;
    double Sn = S + n;
    BBd1 d;
    if (A >= 100.0 || B >= 100.0) {
        double aS = psi_A_rest(n, S);
        double dS = aS + std::log1p(n / S);
        d.lA = (A >= 100.0)
            ? psi_A_rest(y, A) - aS + std::log1p((y * S - A * n) / (A * Sn))
            : psi_diff1(y, A) - dS;
        d.lB = (B >= 100.0)
            ? psi_A_rest(n - y, B) - aS +
                std::log1p(((n - y) * S - B * n) / (B * Sn))
            : psi_diff1(n - y, B) - dS;
        return d;
    }
    double dS = psi_diff1(n, S);
    d.lA = psi_diff1(y, A) - dS;
    d.lB = psi_diff1(n - y, B) - dS;
    return d;
}

inline BBd2 bb_shape_d2(double y, double n, double A, double B) {
    double S = A + B;
    double Sn = S + n;
    BBd2 d;
    double tS  = psi_T_rest(n, S);
    double nSn = n / (S * Sn);
    d.lAA = psi_T_rest(y, A) - tS - y / (A * (A + y)) + nSn;
    d.lBB = psi_T_rest(n - y, B) - tS -
        (n - y) / (B * (B + n - y)) + nSn;
    d.lAB = -tS + nSn;
    return d;
}

// psi'(x + k) - psi'(x): for an integer shift up to 32 minus the sum of
// 1/(x + i)^2, i < k; otherwise two trigammas below the series' crossover of
// psi_T_rest(), the series minus k/(x(x + k)) above it
inline double psi_diff2(double k, double x) {
    if (x >= 100.0) return psi_T_rest(k, x) - k / (x * (x + k));
    if (k <= 32.0 && k >= 0.0 && k == std::floor(k)) {
        double s = 0.0;
        for (int i = 0; i < (int) k; ++i) {
            double r = 1.0 / (x + i);
            s += r * r;
        }
        return -s;
    }
    return R::trigamma(x + k) - R::trigamma(x);
}

// bb_shape_d2() by the same split as bb_shape_d1_mixed(): the form of
// bb_shape_d2() for a component whose shape is at or above 100, the
// differences of psi_diff2() below it
inline BBd2 bb_shape_d2_mixed(double y, double n, double A, double B) {
    if (A >= 100.0 && B >= 100.0) return bb_shape_d2(y, n, A, B);
    double S = A + B;
    double Sn = S + n;
    BBd2 d;
    double dS = psi_diff2(n, S);
    if (A >= 100.0 || B >= 100.0) {
        double tS = psi_T_rest(n, S);
        double nSn = n / (S * Sn);
        d.lAA = (A >= 100.0)
            ? psi_T_rest(y, A) - tS - y / (A * (A + y)) + nSn
            : psi_diff2(y, A) - dS;
        d.lBB = (B >= 100.0)
            ? psi_T_rest(n - y, B) - tS - (n - y) / (B * (B + n - y)) + nSn
            : psi_diff2(n - y, B) - dS;
        d.lAB = -dS;
        return d;
    }
    d.lAA = psi_diff2(y, A) - dS;
    d.lBB = psi_diff2(n - y, B) - dS;
    d.lAB = -dS;
    return d;
}

// both orders, for the components in (mu, sigma), whose second derivatives
// read the first through the curvature of the map
inline BBderiv bb_shape_derivs(double y, double n, double A, double B) {
    const BBd1 d1 = bb_shape_d1(y, n, A, B);
    const BBd2 d2 = bb_shape_d2(y, n, A, B);
    BBderiv d;
    d.lA = d1.lA;
    d.lB = d1.lB;
    d.lAA = d2.lAA;
    d.lBB = d2.lBB;
    d.lAB = d2.lAB;
    return d;
}

// The shapes and their first two derivatives in (mu, sigma).
struct BBmap {
    double A, B;
    double Am, As, Amm, Ams, Ass;
    double Bm, Bs, Bmm, Bms, Bss;
};

inline BBmap bb_map(double mu, double s) {
    BBmap m;
    double s2 = s * s, s3 = s2 * s;
    m.A = mu / s;          m.B = (1.0 - mu) / s;
    m.Am = 1.0 / s;        m.Bm = -1.0 / s;
    m.As = -mu / s2;       m.Bs = -(1.0 - mu) / s2;
    m.Amm = 0.0;           m.Bmm = 0.0;
    m.Ams = -1.0 / s2;     m.Bms = 1.0 / s2;
    m.Ass = 2.0 * mu / s3; m.Bss = 2.0 * (1.0 - mu) / s3;
    return m;
}

const double kBBfk[5] = {0.0, 1.0, -1.0, 2.0, -6.0};
const int kBBpu[3] = {0, 1, 0}, kBBpv[3] = {0, 1, 1};

// the cumulative power sums to order K (3 or 4) and log P(Y = 0)
struct BBSums {
    int N, K;
    std::vector<double> SA, SB;
    double SC[5];
    double lp0;
    BBSums(double a, double b, int N_, int K_)
        : N(N_), K(K_), SA(K_ * (N_ + 1), 0.0), SB(K_ * (N_ + 1), 0.0),
          lp0(0.0) {
        for (int k = 0; k < 5; ++k) SC[k] = 0.0;
        for (int m = 0; m < N; ++m) {
            double ia = 1.0 / (a + m), ib = 1.0 / (b + m), ic = 1.0 / (a + b + m);
            double pa = 1.0, pb = 1.0, pc = 1.0;
            for (int k = 1; k <= K; ++k) {
                pa *= ia; pb *= ib; pc *= ic;
                SA[(k - 1) * (N + 1) + m + 1] = SA[(k - 1) * (N + 1) + m] + pa;
                SB[(k - 1) * (N + 1) + m + 1] = SB[(k - 1) * (N + 1) + m] + pb;
                SC[k] += pc;
            }
            lp0 -= std::log1p(a / (b + m));
        }
    }
    // the derivatives at y = j with na differentiations in a and nb in b,
    // na + nb <= K, into D
    void at(int j, double D[5][5]) const {
        for (int na = 0; na <= K; ++na) {
            for (int nb = 0; na + nb <= K; ++nb) {
                int k = na + nb;
                if (k == 0) { D[na][nb] = 0.0; continue; }
                double c = -kBBfk[k] * SC[k];
                if (nb == 0) c += kBBfk[k] * SA[(k - 1) * (N + 1) + j];
                else if (na == 0) c += kBBfk[k] * SB[(k - 1) * (N + 1) + N - j];
                D[na][nb] = c;
            }
        }
    }
};

inline double bb_next_lp(double lp, int j, int N, double a, double b) {
    return lp + std::log((N - j) / (j + 1.0)) +
        std::log((j + a) / (N - j - 1.0 + b));
}


// ---- betabinom1 (mu, sigma): the components ---------------------------------

// D carries lA and lB: a BBd1 or a BBderiv
template <class D>
inline double betabinom1_score_mu(const D& d, const BBmap& mp) {
    return d.lA * mp.Am + d.lB * mp.Bm;
}

template <class D>
inline double betabinom1_score_sigma(const D& d, const BBmap& mp) {
    return d.lA * mp.As + d.lB * mp.Bs;
}

inline double betabinom1_hess_mu_mu(const BBderiv& d, const BBmap& mp) {
    return d.lAA * mp.Am * mp.Am
         + 2.0 * d.lAB * mp.Am * mp.Bm
         + d.lBB * mp.Bm * mp.Bm
         + d.lA * mp.Amm + d.lB * mp.Bmm;
}

inline double betabinom1_hess_mu_sigma(const BBderiv& d, const BBmap& mp) {
    return d.lAA * mp.Am * mp.As
         + d.lAB * (mp.Am * mp.Bs + mp.As * mp.Bm)
         + d.lBB * mp.Bm * mp.Bs
         + d.lA * mp.Ams + d.lB * mp.Bms;
}

inline double betabinom1_hess_sigma_sigma(const BBderiv& d, const BBmap& mp) {
    return d.lAA * mp.As * mp.As
         + 2.0 * d.lAB * mp.As * mp.Bs
         + d.lBB * mp.Bs * mp.Bs
         + d.lA * mp.Ass + d.lB * mp.Bss;
}

// E[l_kk] by the exact sum over the support, k = 0 for mu and 1 for sigma
inline double betabinom1_expected_kk(int k, double m, double s, double size) {
    BBmap mp = bb_map(m, s);
    int N = (int) size;
    double e = 0.0;
    for (int j = 0; j <= N; j++) {
        double p = std::exp(bb_log_mass(j, mp.A, mp.B, size));
        BBderiv d = bb_shape_derivs(j, size, mp.A, mp.B);
        e += p * (k == 0 ? betabinom1_hess_mu_mu(d, mp)
                         : betabinom1_hess_sigma_sigma(d, mp));
    }
    return e;
}

// ---- the shapes: the expected information and its first derivatives --------

// E[l_ab] for the pairs (aa, bb, ab) in o[0..2] and d_c E[l_ab] in
// o[3 + 2 r + c], r the pair and c the shape differentiated, by one pass
// over the support with BBSums
inline void bb_shapes_dexp1(double a, double b, int N, double* o) {
    const BBSums S(a, b, N, 3);
    double lp = S.lp0;
    double acc[9];
    for (int k = 0; k < 9; ++k) acc[k] = 0.0;
    for (int j = 0; j <= N; ++j) {
        double D[5][5];
        S.at(j, D);
        double p = std::exp(lp);
        int w = 0;
        for (int r = 0; r < 3; ++r) {
            int n0 = (kBBpu[r] == 0) + (kBBpv[r] == 0), n1 = 2 - n0;
            acc[w++] += p * D[n0][n1];
        }
        for (int r = 0; r < 3; ++r) {
            int n0 = (kBBpu[r] == 0) + (kBBpv[r] == 0), n1 = 2 - n0;
            double lab = D[n0][n1];
            for (int c = 0; c < 2; ++c) {
                int c0 = n0 + (c == 0), c1 = n1 + (c == 1);
                double lc = D[c == 0][c == 1];
                acc[w++] += p * (D[c0][c1] + lab * lc);
            }
        }
        if (j < N) lp = bb_next_lp(lp, j, N, a, b);
    }
    for (int k = 0; k < 9; ++k) o[k] = acc[k];
}

// E[l_kk] and d_k E[l_kk] in the shapes, k = 0 for alpha and 1 for beta:
// entries k and 3 + 3k of bb_shapes_dexp1(), each accumulated alone with the
// same terms in the same order
inline void bb_shapes_info_diag(int k, double a, double b, int N,
                                double* out) {
    const BBSums S(a, b, N, 3);
    double lp = S.lp0;
    const int n0 = (k == 0) ? 2 : 0, n1 = 2 - n0;
    const int c0 = n0 + (k == 0), c1 = n1 + (k == 1);
    double e = 0.0, de = 0.0;
    for (int j = 0; j <= N; ++j) {
        double D[5][5];
        S.at(j, D);
        double p = std::exp(lp);
        double lab = D[n0][n1];
        double lc = D[k == 0][k == 1];
        e += p * lab;
        de += p * (D[c0][c1] + lab * lc);
        if (j < N) lp = bb_next_lp(lp, j, N, a, b);
    }
    out[0] = e;
    out[1] = de;
}

// ---- betabinom1: the shapes' derivatives carried to (mu, sigma) -------------

// d_c E[l_ab] in (mu, sigma), parameters indexed 0 = mu and 1 = sigma, from
// the shapes' o[9] of bb_shapes_dexp1() at a = mu/sigma, b = (1 - mu)/sigma,
// by the chain dexpected_chain() writes in R:
//   sum_ij [sum_k dE_ijk J_kc] J_ia J_jb + E_ij (J_iac J_jb + J_ia J_jbc),
// with J the partials of the shapes in (mu, sigma)
inline double betabinom1_dexpected_entry(int pa, int pb, int pc, double m,
                                         double s, const double* o) {
    const double q = 1.0 - m;
    const double r1 = -1.0 / (s * s);
    const double r2 = 2.0 / R_pow(s, 3.0);
    double J1[2][2], J2[2][2][2];
    J1[0][0] = 1.0 / s;    J1[1][0] = -1.0 / s;
    J1[0][1] = m * r1;     J1[1][1] = q * r1;
    J2[0][0][0] = 0.0;     J2[1][0][0] = 0.0;
    J2[0][0][1] = r1;      J2[0][1][0] = r1;
    J2[1][0][1] = -r1;     J2[1][1][0] = -r1;
    J2[0][1][1] = m * r2;  J2[1][1][1] = q * r2;
    const int pr[2][2] = {{0, 2}, {2, 1}};
    double acc = 0.0;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            const int r = pr[i][j];
            const double Eij = o[r];
            const double dEc = o[3 + 2 * r + 0] * J1[0][pc] +
                o[3 + 2 * r + 1] * J1[1][pc];
            acc = acc + dEc * J1[i][pa] * J1[j][pb] +
                Eij * (J2[i][pa][pc] * J1[j][pb] + J1[i][pa] * J2[j][pb][pc]);
        }
    }
    return acc;
}

// ---- the routers of the scalar registry -------------------------------------

// th = (mu, sigma, size)
inline void betabinom1_score_curv(int k, double y, const double* th,
                                  double* out) {
    const double m = th[0], s = th[1], size = th[2];
    const BBmap mp = bb_map(m, s);
    const BBderiv d = bb_shape_derivs(y, size, mp.A, mp.B);
    if (k == 0) {
        out[0] = betabinom1_score_mu(d, mp);
        out[1] = betabinom1_hess_mu_mu(d, mp);
    } else {
        out[0] = betabinom1_score_sigma(d, mp);
        out[1] = betabinom1_hess_sigma_sigma(d, mp);
    }
}

inline void betabinom1_info_dinfo(int k, double y, const double* th,
                                  double* out) {
    (void) y;
    const double m = th[0], s = th[1], size = th[2];
    out[0] = betabinom1_expected_kk(k, m, s, size);
    double o[9];
    bb_shapes_dexp1(m / s, (1.0 - m) / s, (int) size, o);
    out[1] = betabinom1_dexpected_entry(k, k, k, m, s, o);
}

// th = (alpha, beta, size)
inline void betabinom2_score_curv(int k, double y, const double* th,
                                  double* out) {
    const BBd1 d1 = bb_shape_d1_mixed(y, th[2], th[0], th[1]);
    const BBd2 d2 = bb_shape_d2_mixed(y, th[2], th[0], th[1]);
    if (k == 0) {
        out[0] = d1.lA;
        out[1] = d2.lAA;
    } else {
        out[0] = d1.lB;
        out[1] = d2.lBB;
    }
}

inline void betabinom2_info_dinfo(int k, double y, const double* th,
                                  double* out) {
    (void) y;
    bb_shapes_info_diag(k, th[0], th[1], (int) th[2], out);
}

} // namespace d7

#endif
