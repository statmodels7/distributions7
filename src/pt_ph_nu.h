#ifndef D7_PT_PH_NU_H
#define D7_PT_PH_NU_H

#include <cmath>
#include <cstdlib>
#include "pt_bessel_k.h"
#include "pt_sqrt.h"

// The nu-only part of the two pseudo-Huber families: with x = sqrt(nu),
//   B(nu) = log(x K_1(x))  and  R(nu) = x K_2(x)/K_1(x),
// and their derivatives in nu to order five. The pseudo-Huber log-density
// carries -B(nu), the standard-deviation family R(nu) and
// h(nu) = log(R)/2 - B(nu).
//
// Three routes, each where it keeps its digits:
//   nu < 1: the series of K_1 and K_2 at the origin. With u = nu/4,
//     T_n(u) = 2 (x/2)^n K_n(x) = sum_(k<n) (n-k-1)!/k! (-u)^k
//       + (-1)^(n+1) log(u) sum_k u^(n+k)/(k! (n+k)!)
//       + (-1)^n sum_k (psi(k+1) + psi(n+k+1)) u^(n+k)/(k! (n+k)!),
//     so that B = log T_1 and R = 2 T_2/T_1; the derivatives of log T_n
//     come from those of T_n, which the log terms keep large, so no two
//     terms cancel.
//   1 <= nu < 400: K_p(x) = (x/2)^p/2 int t^(-p-1) exp(-t - nu/(4t)) dt, so
//     that the derivatives of log J_p(nu) = log int t^(-p-1) exp(-t - nu/(4t))
//     dt are (-1)^k kappa_k(W), the cumulants of W = 1/(4T) under the
//     generalized inverse Gaussian density of the integrand. They are formed
//     from central moments, by the trapezoidal rule in log t centered at the
//     mode (exponentially convergent for this analytic integrand); the
//     derivatives of log K_n from the ratios K_m/K_n lost up to 1e-10 at
//     order five in this band.
//   nu >= 400: the asymptotic series K_n(x) ~ sqrt(pi/(2x)) e^(-x)
//     sum_k a_k(n) x^(-k). With t = 1/x, B = log(pi/2)/2 + log(nu)/4 - x +
//     F_1(t) and log R = log(nu)/2 + F_2(t) - F_1(t), F_n the power series
//     of log(sum_k a_k(n) t^k); each power t^m = nu^(-m/2) is differentiated
//     exactly, so the leading terms that the ratio route subtracts are never
//     formed.

namespace d7 {

const int PH_ASY_M = 40;

struct PhAsy {
    double f1[PH_ASY_M + 1], g[PH_ASY_M + 1];
    PhAsy() {
        double a1[PH_ASY_M + 1], a2[PH_ASY_M + 1];
        a1[0] = a2[0] = 1.0;
        for (int k = 1; k <= PH_ASY_M; k++) {
            const double o = (2.0 * k - 1.0) * (2.0 * k - 1.0);
            a1[k] = a1[k - 1] * (4.0 - o) / (8.0 * k);
            a2[k] = a2[k - 1] * (16.0 - o) / (8.0 * k);
        }
        double l1[PH_ASY_M + 1], l2[PH_ASY_M + 1];
        logser(a1, l1);
        logser(a2, l2);
        for (int m = 0; m <= PH_ASY_M; m++) {
            f1[m] = l1[m];
            g[m] = l2[m] - l1[m];
        }
    }
    // log of 1 + sum_k a_k t^k as a power series
    static void logser(const double* a, double* l) {
        l[0] = 0.0;
        for (int m = 1; m <= PH_ASY_M; m++) {
            double s = 0.0;
            for (int k = 1; k < m; k++) s += k * l[k] * a[m - k];
            l[m] = a[m] - s / m;
        }
    }
};

inline const PhAsy& ph_asy() {
    static const PhAsy A;
    return A;
}

// d^k/dnu^k of log(nu): (-1)^(k-1) (k-1)!/nu^k
inline double ph_dlog(double v, int k) {
    double f = 1.0;
    for (int i = 2; i < k; i++) f *= i;
    return ((k % 2) ? f : -f) / std::pow(v, k);
}

// c_k = d^k/dnu^k nu^(1/2) times nu^(k - 1/2)
inline double ph_csqrt(int k) {
    double c = 1.0;
    for (int i = 0; i < k; i++) c *= 0.5 - i;
    return c;
}

// the derivatives of log from those of the function: c_k from a_j = F^(j)/F
inline void ph_logder(const double* a, int K, double* c) {
    for (int k = 1; k <= K; k++) {
        double s = a[k], binom = 1.0;
        for (int j = 1; j < k; j++) {
            s -= binom * c[j] * a[k - j];
            binom = binom * (k - j) / j;
        }
        c[k] = s;
    }
}

// T_n^(j)(u), j = 0..K, by the series at the origin
inline void ph_tn_series(int n, double u, int K, double* T) {
    const double lu = std::log(u);
    const double euler = 0.57721566490153286060651209008240243;
    for (int j = 0; j <= K; j++) T[j] = 0.0;
    // the polynomial head
    for (int k = 0; k < n; k++) {
        double c = 1.0;
        for (int i = 2; i <= n - k - 1; i++) c *= i;
        for (int i = 2; i <= k; i++) c /= i;
        if (k % 2) c = -c;
        // d^j u^k
        for (int j = 0; j <= K && j <= k; j++) {
            double ff = 1.0;
            for (int i = 0; i < j; i++) ff *= (k - i);
            T[j] += c * ff * std::pow(u, k - j);
        }
    }
    const double sgn_log = (n % 2) ? 1.0 : -1.0;   // (-1)^(n+1)
    double fk = 1.0, fnk = 1.0;                     // k!, (n+k)!
    for (int i = 2; i <= n; i++) fnk *= i;
    double hk = 0.0, hnk = 0.0;                     // harmonic numbers H_k, H_(n+k)
    for (int i = 1; i <= n; i++) hnk += 1.0 / i;
    for (int k = 0; k < 40; k++) {
        if (k > 0) {
            fk *= k;
            fnk *= (n + k);
            hk += 1.0 / k;
            hnk += 1.0 / (n + k);
        }
        const int m = n + k;
        const double ca = 1.0 / (fk * fnk);
        const double cb = (hk + hnk - 2.0 * euler) * ca;   // psi(k+1) + psi(n+k+1)
        double term_size = 0.0;
        for (int j = 0; j <= K; j++) {
            // d^j (u^m log u) and d^j u^m
            double ffj = 1.0;
            for (int i = 0; i < j; i++) ffj *= (m - i);
            double rest = 0.0, ff = 1.0, binom = 1.0;
            for (int i = 0; i < j; i++) {
                // C(j, i) (m)_i (-1)^(j-i-1) (j-i-1)!
                double fr = 1.0;
                for (int q = 2; q <= j - i - 1; q++) fr *= q;
                rest += binom * ff * (((j - i - 1) % 2) ? -fr : fr);
                ff *= (m - i);
                binom = binom * (j - i) / (i + 1);
            }
            const double up = std::pow(u, m - j);
            const double dlog = up * (ffj * lu + rest);
            const double dpow = ffj * up;
            const double add = sgn_log * ca * dlog - sgn_log * cb * dpow;
            T[j] += add;
            if (j == 0) term_size = std::fabs(add);
        }
        if (k > 2 && term_size < 1e-18 * std::fabs(T[0])) break;
    }
}

// Disabled in T2b (2026-10-11): no caller since the band 1 <= nu < 400
// reads the cumulants of ph_gig_cumulants(); this route lost up to 1e-10 at
// order five there. Kept for reference.
#if 0
// the derivatives of log K_n(sqrt(nu)) in nu, n = 1, 2, by the Bessel ratios
inline void ph_ratio_route(double v, int K, double L[3][6], double& R0) {
    const double t = d7::sqrt_cr(v);
    double km[9] = {0};
    for (int m = 0; m <= K + 2; m++) km[m] = bessel_k_scaled(t, (double) m);
    R0 = t * km[2] / km[1];
    double tk[6] = {0};
    {
        double c = 1.0;
        for (int k = 1; k <= K; k++) {
            c *= (0.5 - (k - 1));
            tk[k] = c * std::pow(v, 0.5 - k);
        }
    }
    for (int n = 1; n <= 2; n++) {
        double rho[6] = {0};
        for (int j = 1; j <= K; j++) {
            double s = 0.0, binom = 1.0;
            for (int i = 0; i <= j; i++) {
                if (i > 0) binom = binom * (j - i + 1) / i;
                s += binom * km[std::abs(n - j + 2 * i)];
            }
            rho[j] = std::pow(-0.5, j) * s / km[n];
        }
        double g[6] = {0};
        const double r1 = rho[1];
        g[1] = r1;
        if (K >= 2) g[2] = rho[2] - r1 * r1;
        if (K >= 3) g[3] = rho[3] - 3.0 * r1 * rho[2] + 2.0 * r1 * r1 * r1;
        if (K >= 4) g[4] = rho[4] - 4.0 * r1 * rho[3] - 3.0 * rho[2] * rho[2] +
                           12.0 * r1 * r1 * rho[2] - 6.0 * std::pow(r1, 4);
        if (K >= 5) g[5] = rho[5] - 5.0 * r1 * rho[4] - 10.0 * rho[2] * rho[3] +
                           20.0 * r1 * r1 * rho[3] + 30.0 * r1 * rho[2] * rho[2] -
                           60.0 * std::pow(r1, 3) * rho[2] + 24.0 * std::pow(r1, 5);
        L[n][1] = g[1] * tk[1];
        if (K >= 2) L[n][2] = g[2] * tk[1] * tk[1] + g[1] * tk[2];
        if (K >= 3) L[n][3] = g[3] * std::pow(tk[1], 3) + 3.0 * g[2] * tk[1] * tk[2] +
                              g[1] * tk[3];
        if (K >= 4) L[n][4] = g[4] * std::pow(tk[1], 4) +
                              6.0 * g[3] * tk[1] * tk[1] * tk[2] +
                              g[2] * (3.0 * tk[2] * tk[2] + 4.0 * tk[1] * tk[3]) +
                              g[1] * tk[4];
        if (K >= 5) L[n][5] = g[5] * std::pow(tk[1], 5) +
                              10.0 * g[4] * std::pow(tk[1], 3) * tk[2] +
                              g[3] * (15.0 * tk[1] * tk[2] * tk[2] +
                                      10.0 * tk[1] * tk[1] * tk[3]) +
                              g[2] * (10.0 * tk[2] * tk[3] + 5.0 * tk[1] * tk[4]) +
                              g[1] * tk[5];
    }
}
#endif

// kap[k], k = 1..K: the cumulants of W = 1/(4T), T with density
// proportional to t^(-p-1) exp(-t - nu/(4t)) on t > 0
inline void ph_gig_cumulants(int p, double v, int K, double* kap) {
    // the mode of the integrand in s = log t, where t^-p exp(-t - nu/(4t))
    // peaks, and its width there
    const double ts = 0.5 * (-p + std::sqrt((double) p * p + v));
    const double sd = 1.0 / std::sqrt(ts + 0.25 * v / ts);
    // the density is skewed at small nu: a finer rule there
    const bool fine = v < 30.0;
    const double s0 = std::log(ts), h = sd / (fine ? 6.0 : 2.5);
    const int N = fine ? 90 : 30;   // nodes on each side, 15 or 12 widths
    const double l0 = -p * s0 - ts - 0.25 * v / ts;
    double w[181], W[181], sw = 0.0, mu = 0.0;
    for (int i = -N; i <= N; i++) {
        const double s = s0 + i * h, t = std::exp(s);
        const double lw = -p * s - t - 0.25 * v / t - l0;
        w[i + N] = std::exp(lw);
        W[i + N] = 0.25 / t;
        sw += w[i + N];
        mu += w[i + N] * W[i + N];
    }
    mu /= sw;
    double m[6] = {0, 0, 0, 0, 0, 0};
    for (int i = 0; i <= 2 * N; i++) {
        const double d = W[i] - mu;
        double dk = d;
        for (int k = 2; k <= K; k++) {
            dk *= d;
            m[k] += w[i] * dk;
        }
    }
    for (int k = 2; k <= K; k++) m[k] /= sw;
    kap[1] = mu;
    if (K >= 2) kap[2] = m[2];
    if (K >= 3) kap[3] = m[3];
    if (K >= 4) kap[4] = m[4] - 3.0 * m[2] * m[2];
    if (K >= 5) kap[5] = m[5] - 10.0 * m[3] * m[2];
}

// dB[k] = B^(k)(nu), dlR[k] = (log R)^(k)(nu), k = 1..K (K <= 5), and R(nu)
inline void ph_nu_derivs_raw(double v, int K, double* dB, double* dlR, double& R0) {
    if (v < 1.0) {
        const double u = 0.25 * v;
        double T1[6], T2[6], a1[6], a2[6], c1[6], c2[6];
        ph_tn_series(1, u, K, T1);
        ph_tn_series(2, u, K, T2);
        for (int j = 1; j <= K; j++) {
            a1[j] = T1[j] / T1[0];
            a2[j] = T2[j] / T2[0];
        }
        ph_logder(a1, K, c1);
        ph_logder(a2, K, c2);
        double sc = 1.0;
        for (int k = 1; k <= K; k++) {
            sc *= 0.25;
            dB[k] = sc * c1[k];
            dlR[k] = sc * (c2[k] - c1[k]);
        }
        R0 = 2.0 * T2[0] / T1[0];
        return;
    }
    if (v >= 400.0) {
        const PhAsy& A = ph_asy();
        const double t = 1.0 / d7::sqrt_cr(v);
        double G0 = 0.0, tm = 1.0;
        for (int m = 0; m <= PH_ASY_M; m++) {
            G0 += A.g[m] * tm;
            tm *= t;
        }
        R0 = d7::sqrt_cr(v) * std::exp(G0);
        for (int k = 1; k <= K; k++) {
            double s1 = 0.0, sg = 0.0;
            tm = 1.0;
            for (int m = 0; m <= PH_ASY_M; m++) {
                double e = 1.0;
                for (int q = 0; q < k; q++) e *= (-0.5 * m - q);
                s1 += A.f1[m] * tm * e;
                sg += A.g[m] * tm * e;
                tm *= t;
            }
            const double vk = std::pow(v, k);
            dB[k] = 0.25 * ph_dlog(v, k) - ph_csqrt(k) * std::pow(v, 0.5 - k) + s1 / vk;
            dlR[k] = 0.5 * ph_dlog(v, k) + sg / vk;
        }
        return;
    }
    // B = log(nu/4) + log J_1 and log R = log(nu/4) + log 2 + log J_2 - log J_1
    {
        const double t = d7::sqrt_cr(v);
        R0 = t * bessel_k_scaled(t, 2.0) / bessel_k_scaled(t, 1.0);
    }
    double k1[6], k2[6];
    ph_gig_cumulants(1, v, K, k1);
    ph_gig_cumulants(2, v, K, k2);
    for (int k = 1; k <= K; k++) {
        const double sg = (k % 2) ? -1.0 : 1.0;
        dB[k] = ph_dlog(v, k) + sg * k1[k];
        dlR[k] = ph_dlog(v, k) + sg * (k2[k] - k1[k]);
    }
}

// the same, with the last nu kept per thread: the kernels and the
// quadrature rules evaluate one nu over runs of consecutive rows, and the
// middle route costs two quadratures
inline void ph_nu_derivs(double v, int K, double* dB, double* dlR, double& R0) {
    struct Last { double v = -1.0; int K = 0; double dB[6], dlR[6], R0; };
    thread_local Last c;
    if (!(v == c.v && K <= c.K)) {
        const int Kc = K > 1 ? K : 1;
        ph_nu_derivs_raw(v, Kc, c.dB, c.dlR, c.R0);
        c.v = v;
        c.K = Kc;
    }
    for (int k = 1; k <= K; k++) {
        dB[k] = c.dB[k];
        dlR[k] = c.dlR[k];
    }
    R0 = c.R0;
}

} // namespace d7

#endif
