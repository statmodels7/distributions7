#ifndef D7_PT_SN2_FDB_H
#define D7_PT_SN2_FDB_H

#include <cmath>
#include "pt_fdb2.h"
#include "pt_logphi.h"

// The G_{a,k} = d_w^a d_g^k F of the centered skew normal by composition,
// for the far left of the step, X < -1, where the generated forms of
// skewnormal2.cpp, written in Z1 = phi/Phi and X, cancel (they lose up to
// 4e-7 at order five). With r = (g/c)^(1/3), Dq^2 = b^2 - (1 - b^2) r^2,
// fA = 1/Dq and fB = (1 + r^2)^(-1/2),
//   F = -log(1 + r^2)/2 - z^2/2 + log Phi(x),
//   z = fB (w + r),  x = A (w + r),  A = r fA fB,
// both linear in w. The derivatives in r of fA, fB and the log are those of
// a power or a logarithm of a quadratic in r, taken by Faa di Bruno; the
// derivatives in g follow from those of r = r(g); and the two outer
// functions, -u^2/2 and log Phi, are composed with z(w, g) and x(w, g) by
// fdb2_compose(), with the derivatives of log Phi from logphi_derivs(),
// which do not cancel below the step.

namespace d7 {

// the incomplete Bell polynomials B_{n,j}(x_1, x_2, ...), n, j <= K
inline void sn2_bell(const double* x, int K, double Bt[6][6]) {
    for (int n = 0; n <= K; n++) for (int j = 0; j <= K; j++) Bt[n][j] = 0.0;
    Bt[0][0] = 1.0;
    for (int n = 1; n <= K; n++) {
        for (int j = 1; j <= n; j++) {
            double s = 0.0, binom = 1.0;   // choose(n - 1, i - 1)
            for (int i = 1; i <= n - j + 1; i++) {
                s += binom * x[i] * Bt[n - i][j - 1];
                binom = binom * (n - i) / i;
            }
            Bt[n][j] = s;
        }
    }
}

// d^n/dr^n of f(u(r)) from f^(j)(u) and the derivatives of u in r
inline void sn2_compose1(const double* fd, const double* ud, int K, double* out) {
    double Bt[6][6];
    sn2_bell(ud, K, Bt);
    out[0] = fd[0];
    for (int n = 1; n <= K; n++) {
        double s = 0.0;
        for (int j = 1; j <= n; j++) s += fd[j] * Bt[n][j];
        out[n] = s;
    }
}

// the Leibniz product of two derivative tables
inline void sn2_leibniz(const double* p, const double* q, int K, double* out) {
    for (int n = 0; n <= K; n++) {
        double s = 0.0, binom = 1.0;
        for (int i = 0; i <= n; i++) {
            s += binom * p[i] * q[n - i];
            binom = binom * (n - i) / (i + 1);
        }
        out[n] = s;
    }
}

// G[a][k] for a + k <= K (K <= 5) at w, at the shape r (r != 0) with Dq
inline void sn2_G_fdb(double w, double r, double Dq, int K, double G[6][6]) {
    const double b2 = 2.0 / M_PI, c = (4.0 - M_PI) / 2.0;
    const double g = c * r * r * r;
    // r(g) and its derivatives in g: r^(q) = r (1/3)_q / g^q
    double rg[6];
    rg[0] = r;
    {
        double ff = 1.0;
        for (int q = 1; q <= K; q++) {
            ff *= (1.0 / 3.0 - (q - 1));
            rg[q] = r * ff / std::pow(g, q);
        }
    }
    // u = Dq^2 = b^2 - (1 - b^2) r^2 and v = 1 + r^2, with their r-derivatives
    double ud[6] = {0, -2.0 * (1.0 - b2) * r, -2.0 * (1.0 - b2), 0, 0, 0};
    double vd[6] = {0, 2.0 * r, 2.0, 0, 0, 0};
    const double v = 1.0 + r * r;
    // f^(j)(u) of u^(-1/2), read with u = Dq^2, and of v^(-1/2), -log(v)/2
    double fa[6], fb[6], fl[6];
    {
        double pf = 1.0;
        for (int j = 0; j <= K; j++) {
            if (j > 0) pf *= (-0.5 - (j - 1));
            fa[j] = pf * std::pow(Dq, -1.0 - 2.0 * j);
            fb[j] = pf * std::pow(v, -0.5 - j);
        }
        fl[0] = -0.5 * std::log(v);
        double f = 1.0;
        for (int j = 1; j <= K; j++) {
            if (j > 1) f *= -(j - 1);
            fl[j] = -0.5 * f / std::pow(v, j);
        }
    }
    double FA[6], FB[6], PL[6];
    sn2_compose1(fa, ud, K, FA);
    sn2_compose1(fb, vd, K, FB);
    sn2_compose1(fl, vd, K, PL);
    // r, r^2 as tables in r, and the products A = r FA FB, Q = r A, Z = r FB
    double R1[6] = {r, 1, 0, 0, 0, 0}, AB[6], A[6], Q[6], Z[6];
    sn2_leibniz(FA, FB, K, AB);
    sn2_leibniz(R1, AB, K, A);
    sn2_leibniz(R1, A, K, Q);
    sn2_leibniz(R1, FB, K, Z);
    // to derivatives in g
    double Bt[6][6];
    sn2_bell(rg, K, Bt);
    auto tog = [&](const double* h, double* out) {
        out[0] = h[0];
        for (int q = 1; q <= K; q++) {
            double s = 0.0;
            for (int j = 1; j <= q; j++) s += h[j] * Bt[q][j];
            out[q] = s;
        }
    };
    double Pg[6], Bg[6], Zg[6], Ag[6], Qg[6];
    tog(PL, Pg);
    tog(FB, Bg);
    tog(Z, Zg);
    tog(A, Ag);
    tog(Q, Qg);
    // the inner tables of z = fB w + r fB and x = A w + Q
    double dz[6][6] = {{0}}, dx[6][6] = {{0}};
    for (int q = 0; q <= K; q++) {
        dz[0][q] = Bg[q] * w + Zg[q];
        dx[0][q] = Ag[q] * w + Qg[q];
        dz[1][q] = Bg[q];
        dx[1][q] = Ag[q];
    }
    // outer derivatives: -u^2/2 at z, log Phi at x
    double fz[6] = {0, -dz[0][0], -1.0, 0, 0, 0};
    double fx[6];
    logphi_derivs(dx[0][0], K, fx);
    for (int a = 0; a <= K; a++) {
        for (int k = 0; a + k <= K; k++) {
            if (a + k == 0) { G[0][0] = 0.0; continue; }
            G[a][k] = (a == 0 ? Pg[k] : 0.0) + fdb2_compose(a, k, fz, dz) +
                fdb2_compose(a, k, fx, dx);
        }
    }
}

} // namespace d7

#endif
