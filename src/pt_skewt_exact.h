#ifndef D7_PT_SKEWT_EXACT_H
#define D7_PT_SKEWT_EXACT_H

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include "pt_sqrt.h"
#include "pt_skewt_table.h"

// The skew t's derivatives in (mu, sigma, alpha, nu) to order five, exact in
// nu. With z = (y - mu)/sigma, rs = sqrt(nu + z^2), a = z/rs, m = nu + 1 and
// w = alpha sqrt(m) a,
//   l = log 2 - log sigma + C(nu) - (nu + 1)/2 log1p(z^2/nu) + log T_m(w),
// C(m) = lgamma((m+1)/2) - lgamma(m/2) - log(m pi)/2. The components are the
// generated closed forms of pt_skewt_table.h in the bounded variables
// (a, 1/rs, 1/sigma, alpha, nu, sqrt(m)), fed with the derivatives of C and
// with B_ij = d_w^i d_m^j log T_m(w).
//
// The B_ij are polynomials in q = t_m(w)/T_m(w), Delta = d_w log t_m(w) - q,
// the partials E_ab of log t_m(w) and the ratios R_0l = (d_m^l T_m(w))/T_m(w)
// (the Riccati form of pt_skewt_table.h). The E_ab are closed forms. The
// R_0l have none: d_m^l T_m(w) is the integral of d_m^l t_m(u) over u, taken
// here by quadrature (skewt_R0()).
//
// The derivatives of C are 2^-k D^(k)(m/2), D(x) = lgamma(x + 1/2) -
// lgamma(x) - log(x)/2 (skewt_cder()). Written as a difference of polygammas
// they lose a factor m to cancellation (2.5e-11 at m = 201); D is shifted to
// x >= 20 by D(x) = D(x + 1) + log((x+1)/x)/2 - log((x + 1/2)/x), whose
// derivatives are expm1/log1p forms, and summed there by its asymptotic
// series, to 7e-16 at every m measured against mpmath.
//
// Each order has its own kernel (skewt_obs<K>()), which computes the
// integrals and the B_ij to its own order and nothing above it.

namespace d7 {

typedef double (*N7Pt)(double, double, int, int);
N7Pt skewt_pt();

// D^(k)(x), k = 1..kmax, into D[1..kmax]
inline void skewt_dder(double x, int kmax, double* D) {
    // c_n of D(x) ~ sum_{n odd} c_n x^-n, from the Bernoulli polynomials:
    // c_n = (-1)^(n+1) (B_{n+1}(1/2) - B_{n+1}(0))/(n (n+1))
    static const double c[10] = {
        -0.125, 0.0052083333333333333333, -0.0015625,
        0.0011858258928571428571, -0.0016818576388888888889,
        0.0038341175426136363636, -0.012819730318509615385,
        0.059100405375162760417, -0.35928737415986902574,
        2.7848617779581170333};
    const int N = (x < 20.0) ? (int) std::ceil(20.0 - x) : 0;
    const double X = x + N, iX = 1.0 / X;
    double sh[6] = {0, 0, 0, 0, 0, 0};
    for (int i = 0; i < N; ++i) {
        const double iy = 1.0 / (x + i);
        const double l1 = std::log1p(iy), lh = std::log1p(0.5 * iy);
        double yk = 1.0;
        for (int k = 1; k <= kmax; ++k) {
            yk *= iy;
            sh[k] += yk * (0.5 * std::expm1(-k * l1) - std::expm1(-k * lh));
        }
    }
    double fact = 1.0;   // (k - 1)!
    for (int k = 1; k <= kmax; ++k) {
        if (k > 1) fact *= (k - 1);
        double as = 0.0;
        for (int t = 9; t >= 0; --t) {
            const int n = 2 * t + 1;
            double pre = 1.0;
            for (int i = 0; i < k; ++i) pre *= (n + i);
            as += c[t] * pre * R_pow_di(iX, n + k);
        }
        if (k & 1) as = -as;
        const double sg = (k & 1) ? 1.0 : -1.0;   // (-1)^(k-1)
        D[k] = as + sg * fact * sh[k];
    }
}

// C^(k)(m) = 2^-k D^(k)(m/2), k = 1..kmax, into C[1..kmax]
inline void skewt_cder(double m, int kmax, double* C) {
    skewt_dder(0.5 * m, kmax, C);
    for (int k = 1; k <= kmax; ++k) C[k] = std::ldexp(C[k], -k);
}

// The m-derivatives of e(u) = log t_m(u), k = 1..KM, into ell[1..KM], with
// Cm[k] = C^(k)(m) and g0 = log1p(u^2/m):
//   ell_k = C^(k)(m) - (m + 1)/2 g_k - k/2 g_(k-1),
//   g_j = d_m^j log1p(u^2/m) = (-1)^(j-1) (j-1)! m^-j (rho^j - 1),
// rho = m/(m + u^2), with rho^j - 1 = -v (1 + rho + ... + rho^(j-1)) and
// v = u^2/(m + u^2), so that no g_j is a difference of near numbers.
template <int KM>
inline void skewt_ell(double u2, double m, double im, double hm1,
                      double g0, const double* Cm, double* ell) {
    const double den = 1.0 / (m + u2), v = u2 * den, rho = m * den;
    double g[KM + 1];
    g[0] = g0;
    double geo = 0.0, rp = 1.0, mj = 1.0, f = 1.0;   // f = (-1)^(j-1) (j-1)!
    for (int j = 1; j <= KM; ++j) {
        geo += rp;
        rp *= rho;
        mj *= im;
        if (j > 1) f *= -(j - 1);
        g[j] = -f * mj * v * geo;
    }
    for (int k = 1; k <= KM; ++k) ell[k] = Cm[k] - (hm1 * g[k] + 0.5 * k * g[k - 1]);
}

// the complete Bell polynomials of ell[1..KM], P[j] = (d_m^j t)/t
template <int KM>
inline void skewt_bell(const double* l, double* P) {
    const double l1 = l[1];
    P[1] = l1;
    if (KM >= 2) P[2] = l[2] + l1 * l1;
    if (KM >= 3) P[3] = l[3] + 3.0 * l1 * l[2] + l1 * l1 * l1;
    if (KM >= 4) {
        const double l12 = l1 * l1;
        P[4] = l[4] + 4.0 * l1 * l[3] + 3.0 * l[2] * l[2] + 6.0 * l12 * l[2] +
            l12 * l12;
    }
    if (KM >= 5) {
        const double l12 = l1 * l1;
        P[5] = l[5] + 5.0 * l1 * l[4] + 10.0 * l[2] * l[3] + 10.0 * l12 * l[3] +
            15.0 * l1 * l[2] * l[2] + 10.0 * l12 * l1 * l[2] + l12 * l12 * l1;
    }
}

// 12-point Gauss-Legendre on [-1, 1]
const double kSkewtGLx[12] = {
    -0.9815606342467192506905, -0.9041172563704748566785,
    -0.7699026741943046870369, -0.5873179542866174472967,
    -0.3678314989981801937527, -0.1252334085114689154724,
    0.1252334085114689154724, 0.3678314989981801937527,
    0.5873179542866174472967, 0.7699026741943046870369,
    0.9041172563704748566785, 0.9815606342467192506905};
const double kSkewtGLw[12] = {
    0.04717533638651182719462, 0.1069393259953184309603,
    0.1600783285433462263347, 0.2031674267230659217491,
    0.2334925365383548087608, 0.2491470458134027850006,
    0.2491470458134027850006, 0.2334925365383548087608,
    0.2031674267230659217491, 0.1600783285433462263347,
    0.1069393259953184309603, 0.04717533638651182719462};

// 20-point Gauss-Legendre on [-1, 1], for the interval [0, w]
const double kSkewtGL20x[20] = {
    -0.9931285991850949247861, -0.9639719272779137912677,
    -0.9122344282513259058678, -0.8391169718222188233945,
    -0.7463319064601507926143, -0.6360536807265150254528,
    -0.5108670019508270980044, -0.3737060887154195606725,
    -0.2277858511416450780805, -0.07652652113349733375464,
    0.07652652113349733375464, 0.2277858511416450780805,
    0.3737060887154195606725, 0.5108670019508270980044,
    0.6360536807265150254528, 0.7463319064601507926143,
    0.8391169718222188233945, 0.9122344282513259058678,
    0.9639719272779137912677, 0.9931285991850949247861};
const double kSkewtGL20w[20] = {
    0.01761400713915211831186, 0.04060142980038694133104,
    0.06267204833410906356951, 0.08327674157670474872476,
    0.1019301198172404350368, 0.1181945319615184173124,
    0.1316886384491766268985, 0.1420961093183820513293,
    0.1491729864726037467878, 0.1527533871307258506981,
    0.1527533871307258506981, 0.1491729864726037467878,
    0.1420961093183820513293, 0.1316886384491766268985,
    0.1181945319615184173124, 0.1019301198172404350368,
    0.08327674157670474872476, 0.06267204833410906356951,
    0.04060142980038694133104, 0.01761400713915211831186};

// The tail's nodes in tau, their weights and the factors exp(tau/s) at them
// depend on s = min(m, 8) alone, so they are computed once per s and kept
// per thread: 96 exponentials an observation fewer.
struct SkewtTailNodes {
    double sc = -1.0;
    double hw[96], es[96];   // weight times panel half-width; exp(tau/s)/s
    double em1[96];          // expm1(tau/s)
    void at(double s) {
        if (s == sc) return;
        sc = s;
        static const double br[9] = {0, 0.5, 1, 2, 4, 8, 16, 32, 64};
        for (int p = 0; p < 8; ++p) {
            const double lo = br[p], hh = 0.5 * (br[p + 1] - lo), mid = lo + hh;
            for (int i = 0; i < 12; ++i) {
                const int k = 12 * p + i;
                const double tau = mid + hh * kSkewtGLx[i];
                hw[k] = kSkewtGLw[i] * hh;
                em1[k] = std::expm1(tau / s);
                es[k] = (em1[k] + 1.0) / s;
            }
        }
    }
};

// R0[j] = delta_j = (d_m^j T_m(w))/T_m(w) - P_j(w), j = 1..KM, with
// q = t_m(w)/T_m(w) and Pw[j] = P_j(w) = (d_m^j t_m(w))/t_m(w).
//
// T_m(0) = 1/2 at every m, so d_m^j T_m(w) = int_0^w d_m^j t_m(u) du, a
// short interval for |w| <= 2, taken by 20-point Gauss-Legendre, where the
// integral over the tail would be a difference of two halves of a vanishing
// total (3e-9 at m = 31, w = 0.3). Beyond, it is minus the integral over the
// tail on the far side of w, in u = w -/+ |w| expm1(tau/s), s = min(m, 8),
// where the density's algebraic tail decays like exp(-tau): 12-point
// Gauss-Legendre on the panels 0, 0.5, 1, 2, ..., 64 in tau, stopped once a
// panel adds less than 1e-20 of the tail's mass, which does not depend on
// KM, so that the sum for one j is the same number at every order. On the
// left the integrand is t(u) (P_j(u) - P_j(w)), whose integral is T delta_j
// directly, the tail's mass being T: R_0j and P_j(w) nearly coincide in the
// left tail, and their difference formed after the integral lost 1e3 at
// w = -78. Against mpmath at 150 digits, over m from 1.3 to 201 and w from
// -1000 to 1000, every ratio is within 2.6e-15 for m <= 5, 2.2e-14 at m = 31
// and 1.1e-12 at m = 201.
template <int KM>
inline void skewt_R0(double w, double m, double q, const double* Cm,
                     const double* Pw, double* R0) {
    const double im = 1.0 / m, hm1 = 0.5 * (m + 1.0);
    const double gw = std::log1p(w * w * im);
    double acc[KM + 1];
    for (int j = 1; j <= KM; ++j) acc[j] = 0.0;
    double ell[KM + 1], P[KM + 1];
    if (std::fabs(w) <= 2.0) {
        const double hw = 0.5 * w;
        for (int i = 0; i < 20; ++i) {
            const double u = hw * (kSkewtGL20x[i] + 1.0), u2 = u * u;
            const double g0 = std::log1p(u2 * im);
            const double r = kSkewtGL20w[i] * std::exp(-hm1 * (g0 - gw));
            skewt_ell<KM>(u2, m, im, hm1, g0, Cm, ell);
            skewt_bell<KM>(ell, P);
            for (int j = 1; j <= KM; ++j) acc[j] += r * P[j];
        }
        for (int j = 1; j <= KM; ++j) R0[j] = q * hw * acc[j] - Pw[j];
        return;
    }
    thread_local SkewtTailNodes nodes;
    nodes.at(std::min(m, 8.0));
    const double side = (w < 0.0) ? -1.0 : 1.0;
    const double c0 = std::fabs(w), sc0 = side * c0;
    double mass = 0.0;
    for (int p = 0; p < 8; ++p) {
        double add[KM + 1], add0 = 0.0;
        for (int j = 1; j <= KM; ++j) add[j] = 0.0;
        for (int i = 0; i < 12; ++i) {
            const int k = 12 * p + i;
            const double u = w + sc0 * nodes.em1[k], u2 = u * u;
            const double g0 = std::log1p(u2 * im);
            const double r = std::exp(-hm1 * (g0 - gw));
            if (!(r > 0.0) || !std::isfinite(u2)) continue;
            const double wt = nodes.hw[k] * c0 * nodes.es[k] * r;
            add0 += wt;
            skewt_ell<KM>(u2, m, im, hm1, g0, Cm, ell);
            skewt_bell<KM>(ell, P);
            if (side < 0.0) for (int j = 1; j <= KM; ++j) P[j] -= Pw[j];
            for (int j = 1; j <= KM; ++j) add[j] += wt * P[j];
        }
        mass += add0;
        for (int j = 1; j <= KM; ++j) acc[j] += add[j];
        if (p >= 3 && add0 <= 1e-20 * mass) break;
    }
    if (side < 0.0) {
        for (int j = 1; j <= KM; ++j) R0[j] = q * acc[j];
    } else {
        for (int j = 1; j <= KM; ++j) R0[j] = -q * acc[j] - Pw[j];
    }
}

// What an observation reads, independent of the order: the bounded
// variables and the quantities of w.
struct SkewtObs {
    double a, irs, isg, rm, m, w, L, q;
};

inline SkewtObs skewt_obs_vars(double y, double mu, double sigma, double al,
                               double nu) {
    SkewtObs o;
    const double z = (y - mu) / sigma;
    if (std::fabs(z) > 1e100) {
        // nu + z^2 overflows near |z| = 1e154; rs = |z| sqrt(1 + nu/z^2)
        const double t = nu / (z * z);
        const double rs = std::fabs(z) * d7::sqrt_cr(1.0 + t);
        o.irs = 1.0 / rs;
        o.a = (z > 0 ? 1.0 : -1.0) / d7::sqrt_cr(1.0 + t);
        o.L = 2.0 * std::log(std::fabs(z)) + std::log1p(t) - std::log(nu);
    } else {
        const double rs = d7::sqrt_cr(nu + z * z);
        o.irs = 1.0 / rs;
        o.a = z * o.irs;
        o.L = std::log1p(z * z / nu);
    }
    o.isg = 1.0 / sigma;
    o.m = nu + 1.0;
    o.rm = d7::sqrt_cr(o.m);
    o.w = al * o.rm * o.a;
    o.q = std::exp(R::dt(o.w, o.m, 1) - skewt_pt()(o.w, o.m, 1, 1));
    return o;
}

// B_ij, i + j <= K, at one observation; Cm[k] = C^(k)(m)
// with_nu = false leaves out the integrals: the components free of nu read
// no B_ij with j >= 1, and are then the same numbers at a fraction of the cost
template <int K>
inline void skewt_btab(const SkewtObs& o, const double* Cm, double B[6][6],
                       bool with_nu = true) {
    double E[6][6] = {{0}};
    double Dl = 0.0;
    // E_0b = ell_b(w), b <= K, and P_0b(w), their complete Bell polynomials
    double ell[K + 1], Pw[K + 1];
    const double w2 = o.w * o.w;
    skewt_ell<K>(w2, o.m, 1.0 / o.m, 0.5 * (o.m + 1.0), std::log1p(w2 / o.m),
                 Cm, ell);
    skewt_bell<K>(ell, Pw);
    for (int b = 1; b <= K; ++b) E[0][b] = ell[b];
    if (K >= 2) {
        // the E_ab with a >= 1 in closed form
        const double den = 1.0 / (o.m + w2);
        // each order fills its own entries, so that an entry is the same
        // number whatever K it is computed for
        skewt_E_2(o.w, den, E);
        if (K >= 3) skewt_E_3(o.w, den, E);
        if (K >= 4) skewt_E_4(o.w, den, E);
        if (K >= 5) skewt_E_5(o.w, den, E);
        // Delta = E_10 - q, the one difference the Riccati form takes
        Dl = E[1][0] - o.q;
    }
    // delta_l = R_0l - P_0l(w), from the integrals
    double G0[6] = {0, 0, 0, 0, 0, 0};
    if (with_nu) skewt_R0<K>(o.w, o.m, o.q, Cm, Pw, G0);
    skewt_B_1(o.q, Dl, E, G0, B);
    if (K >= 2) skewt_B_2(o.q, Dl, E, G0, B);
    if (K >= 3) skewt_B_3(o.q, Dl, E, G0, B);
    if (K >= 4) skewt_B_4(o.q, Dl, E, G0, B);
    if (K >= 5) skewt_B_5(o.q, Dl, E, G0, B);
}

// The derivative constants of one nu: C^(k)(nu) for the density of z and
// C^(k)(nu + 1) for the integrals, k = 1..K
// and computed again only when nu or the order changes
struct SkewtConst {
    double nu = -1.0;
    int K = 0;
    double Cnu[7], Cm[7];
    void at(double v, int k) {
        if (v == nu && k <= K) return;
        nu = v;
        K = k;
        skewt_cder(v, k, Cnu);
        skewt_cder(v + 1.0, k, Cm);
    }
};

// All components of order K at one observation, in the order of
// pt_skewt_table.h
template <int K>
inline void skewt_obs(double y, double mu, double sigma, double al, double nu,
                      const SkewtConst& cc, double* out, bool with_nu = true) {
    const SkewtObs o = skewt_obs_vars(y, mu, sigma, al, nu);
    double B[6][6] = {{0}};
    skewt_btab<K>(o, cc.Cm, B, with_nu);
    switch (K) {
    case 1: skewt_d_1(o.a, o.irs, o.isg, al, nu, o.rm, o.L, cc.Cnu, B, out); break;
    case 2: skewt_d_2(o.a, o.irs, o.isg, al, nu, o.rm, o.L, cc.Cnu, B, out); break;
    case 3: skewt_d_3(o.a, o.irs, o.isg, al, nu, o.rm, o.L, cc.Cnu, B, out); break;
    case 4: skewt_d_4(o.a, o.irs, o.isg, al, nu, o.rm, o.L, cc.Cnu, B, out); break;
    default: skewt_d_5(o.a, o.irs, o.isg, al, nu, o.rm, o.L, cc.Cnu, B, out); break;
    }
}

} // namespace d7

#endif
