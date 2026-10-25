// The per-observation rows of the Poisson-inverse Gaussian's log-mass
// derivatives: pig1_row() and pig2_row() (orders 0 to 4, from pig_hd.cpp) and
// pig1_d5_row() and pig2_d5_row() (order 5, from d5_pig1.cpp and
// d5_pig2.cpp). They live here so that the kernels of pig_hd.cpp and
// d5_pig*.cpp and those of the third and fourth derivatives of the expected
// information (d34e_pig*.cpp) compile the same code.
#ifndef DISTRIBUTIONS7_PT_PIG_ROWS_H
#define DISTRIBUTIONS7_PT_PIG_ROWS_H

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <vector>
#include "d7_par.h"
#include "pig_asym.h"
#include "pt_pig.h"
#include "pt_pig_d5.h"
#include "pt_sqrt.h"

namespace d7pigrow {
using d7::psi_derivs;
using d7::PIG_OFF;
using d7::PIG_LEN;
using d7::pix;
using d7::pig_compose;
using d7::pig_logS_w;
using d7::PigSeq;
using d7::PigSeqW;
using d7::PIG_ASYM_CUT;

inline void pig1_row(double yy, double m, double sg,
                      int order, bool only, double* v,
                      const double* pre = nullptr) {
    // The support test lives here, as it does in every other
    // compiled discrete family: a response that is negative,
    // fractional or not finite has no mass, and the kernel says so
    // rather than leaving a caller to mask afterwards. The value
    // reads -Inf, which is the log-mass there; a derivative reads
    // NaN, there being no derivative to report.
    if (ISNAN(yy)) {
        // a missing response propagates, as it does through dpois()
        int k = only ? PIG_LEN[order] : PIG_OFF[order] + PIG_LEN[order];
        for (int j = 0; j < k; ++j) v[j] = yy;
        return;
    }
    if (!R_finite(yy) || yy < 0.0 || yy != std::floor(yy)) {
        int lo = only ? PIG_OFF[order] : 0;
        int k = only ? PIG_LEN[order] : PIG_OFF[order] + PIG_LEN[order];
        for (int j = 0; j < k; ++j)
            v[j] = (lo + j == 0) ? R_NegInf : R_NaN;
        return;
    }
    double c = 1.0 + 2.0 * sg * m;
    double s = d7::sqrt_cr(c);
    double w = 0.5 * sg / s;
    double L[5];
    if (pre) {
        for (int j = 0; j < 5; ++j) L[j] = pre[j];
    } else {
        pig_logS_w(yy, w, order, L);
    }
    // c is bilinear: c_mu = 2 sigma, c_sigma = 2 mu, c_{mu sigma} = 2
    double C[15] = {0};
    C[0] = c; C[1] = 2.0 * sg; C[2] = 2.0 * m; C[4] = 2.0;
    double ic = 1.0 / c;
    double hlog[5] = {std::log(c), ic, -ic * ic, 2.0 * ic * ic * ic,
                      -6.0 * ic * ic * ic * ic};
    double hsq[5] = {s, 0.5 / s, -0.25 / (s * c), 0.375 / (s * c * c),
                     -0.9375 / (s * c * c * c)};
    double LC[15], S[15];
    pig_compose(hlog, C, order, LC);
    pig_compose(hsq, C, order, S);
    // G = -2 mu phi(s), phi(s) = 1/(1 + s)
    double ip = 1.0 / (1.0 + s);
    double hphi[5] = {ip, -ip * ip, 2.0 * ip * ip * ip,
                      -6.0 * ip * ip * ip * ip, 24.0 * ip * ip * ip * ip * ip};
    // w = (sigma/2) psi(s), psi(s) = 1/s
    double is = 1.0 / s;
    double hpsi[5] = {is, -is * is, 2.0 * is * is * is,
                      -6.0 * is * is * is * is, 24.0 * is * is * is * is * is};
    double PH[15], PS[15], Wp[15], LS[15];
    pig_compose(hphi, S, order, PH);
    pig_compose(hpsi, S, order, PS);
    // Leibniz against the linear factors mu (for G) and sigma (for w)
    double G[15];
    for (int o = 0; o <= order; ++o) {
        for (int j = 0; j <= o; ++j) {
            int i = o - j, k = pix(i, j);
            G[k] = -2.0 * (m * PH[k] + (i > 0 ? i * PH[pix(i - 1, j)] : 0.0));
            Wp[k] = 0.5 * (sg * PS[k] + (j > 0 ? j * PS[pix(i, j - 1)] : 0.0));
        }
    }
    Wp[0] = w;
    pig_compose(L, Wp, order, LS);
    double h = -yy / 2.0;
    double im = 1.0 / m;
    double ym[5] = {yy * std::log(m), yy * im, -yy * im * im,
                    2.0 * yy * im * im * im, -6.0 * yy * im * im * im * im};
    int lo = only ? PIG_OFF[order] : 0;
    int hi = PIG_OFF[order] + PIG_LEN[order];
    for (int o = 0; o <= order; ++o) {
        for (int j = 0; j <= o; ++j) {
            int i = o - j, k = pix(i, j);
            if (k < lo || k >= hi) continue;
            double val = h * LC[k] + G[k] + LS[k];
            if (j == 0) val += ym[i];
            if (k == 0) {
                val -= R::lgammafn(yy + 1.0);
                // the log mass at large counts without the cancellation
                if (!pre && d7::pig_logmass_big_ok(yy, m)) {
                    const double vb = d7::pig1_logmass_big(yy, m, sg);
                    if (std::isfinite(vb)) val = vb;
                }
            }
            v[k - lo] = val;
        }
    }
}

inline void pig2_row(double yy, double m, double aa,
                      int order, bool only, double* v,
                      const double* pre = nullptr) {
    // The support test lives here, as it does in every other
    // compiled discrete family: a response that is negative,
    // fractional or not finite has no mass, and the kernel says so
    // rather than leaving a caller to mask afterwards. The value
    // reads -Inf, which is the log-mass there; a derivative reads
    // NaN, there being no derivative to report.
    if (ISNAN(yy)) {
        // a missing response propagates, as it does through dpois()
        int k = only ? PIG_LEN[order] : PIG_OFF[order] + PIG_LEN[order];
        for (int j = 0; j < k; ++j) v[j] = yy;
        return;
    }
    if (!R_finite(yy) || yy < 0.0 || yy != std::floor(yy)) {
        int lo = only ? PIG_OFF[order] : 0;
        int k = only ? PIG_LEN[order] : PIG_OFF[order] + PIG_LEN[order];
        for (int j = 0; j < k; ++j)
            v[j] = (lo + j == 0) ? R_NegInf : R_NaN;
        return;
    }
        // r = sqrt(m^2 + a^2) and its partials, written in r_m and r_a --
        // which lie in [0, 1] -- against powers of 1/r. The direct form
        // divides by r^3, r^5 and r^7 and multiplies by a^3 and a^4, and
        // a^3 passes double.xmax at alpha = 5.6e102: there r_mma read
        // Inf - Inf and the whole Hessian came back NaN. std::hypot keeps
        // r itself finite where m^2 + a^2 does not.
        double r = std::hypot(m, aa);
        double ir = 1.0 / r, ir2 = ir * ir, ir3 = ir2 * ir;
        // sigma = (m + r) * a^{-2} by the Leibniz rule; u = m + r. The
        // powers of a are formed from ia so that a^4 cannot overflow;
        // where the true partial is below the smallest double the ladder
        // underflows to zero, which is the honest reading there.
        double ia = 1.0 / aa;
        double v0 = ia * ia;
        double v1 = -2.0 * v0 * ia, v2 = 6.0 * v0 * v0,
            v3 = -24.0 * v0 * v0 * ia,
            v4 = 120.0 * v0 * v0 * v0;
        // b = 1/sigma = a / (t + sqrt(1 + t^2)) with t = m/a, which is
        // the same number as a^2/(m + r) with neither a^2 nor m + r formed.
        double t = m * ia, wt = std::hypot(1.0, t);
        double b = aa / (t + wt);
        // F(sigma) = -y log sigma + 1/sigma has F_k = d^k F / dsigma^k
        // carrying b^k, and every Faa di Bruno term of order k has exactly
        // k factors S: the powers cancel between them. The expansion is
        // therefore written in Q_k = F_k / b^k, each linear in b, against
        // T_x = b S_x, each of size 1/a. Formed separately, F_4 = 6 y b^4
        // + 24 b^5 leaves the doubles at alpha = 4.5e61 while the product
        // it belongs to is of size a^{-3}.
        double Q1 = -yy - b;
        double Q2 = yy + 2.0 * b;
        double Q3 = -2.0 * yy - 6.0 * b;
        double Q4 = 6.0 * yy + 24.0 * b;
        double p[5], logS, A1;
        // `pre` carries what psi_derivs() would return, from the recurrence
        // the expected-information kernel runs over the support
        if (pre) {
            for (int j = 0; j < 5; ++j) p[j] = pre[j];
            logS = pre[5];
            A1 = pre[6];
        } else {
            psi_derivs(yy, aa, order, p, &logS, &A1);
        }
        // pure pieces: y log m, -y log(a sigma) and 1/sigma - alpha.
        //
        // a sigma = (m + r)/a = t + sqrt(1 + t^2), so -y log(a sigma) is
        // -y asinh(t): written as -y log a - y log sigma it is a
        // difference of two quantities of size log a.
        //
        // 1/sigma - alpha = a(a - m - r)/(m + r), and alpha - r =
        // -m^2/(alpha + r), so the pair is
        //     -m [a/(m + r)] [1 + m/(a + r)],
        // both factors bounded, tending to -m as alpha grows, which is
        // the Poisson term. Written directly it is a difference of two
        // quantities of size alpha whose value is of size mu: past
        // alpha = 1e16 the log-mass read zero, that is a probability of
        // one, and at 1e18 it read +128, one ulp of 1e18, so the mass
        // over the support summed to Inf. A fit maximizing that reward
        // drove the dispersion out of range instead of stopping.
        double im = 1.0 / m;
        double core = -m * (aa / (m + r)) * (1.0 + m / (aa + r));
        int w = 0;
        v[w + 0] = yy * std::log(m) - yy * std::asinh(t) + core + logS -
            R::lgammafn(yy + 1.0);
        // the log mass at large counts without the cancellation
        if (!pre && d7::pig_logmass_big_ok(yy, m)) {
            const double vb = d7::pig2_logmass_big(yy, m, aa);
            if (std::isfinite(vb)) v[w + 0] = vb;
        }
        if (order == 0) return;
        double r_m = m * ir, r_a = aa * ir;
        double rm2 = r_m * r_m, ra2 = r_a * r_a;
        double u = m + r, u_m = 1.0 + r_m;
        double S_m = u_m * v0;
        double S_a = r_a * v0 + u * v1;
        double T_m = b * S_m, T_a = b * S_a;
        // F(sigma(m, a)) by the same written-out Faa di Bruno
        double P_m = Q1 * T_m;
        // The score in alpha, in closed form for the same reason the value
        // is. Since r^2 - m^2 = alpha^2 the pair is 1/sigma - alpha =
        // (r - m) - alpha, and r - alpha = m^2/(r + alpha), so its
        // derivative is -m^2/(r(r + alpha)) = -r_m m/(r + alpha); the
        // other piece, -y asinh(m/alpha), gives y m/(alpha r). Both are of
        // size alpha^-2, which is what the score is. Assembled as
        // -y/alpha + psi'(alpha) + F_1 S_a the three terms are each of
        // size one: measured, the score was noise past alpha = 1e8 and
        // read exactly -1 past 1e162, where 1/alpha^2 leaves the doubles,
        // so a search saw a slope where the surface is flat.
        double dD_a = -r_m * (m / (r + aa));
        w = only ? 0 : 1;
        v[w + 0] = yy * im + P_m;
        v[w + 1] = yy * m * ir * ia + A1 + dD_a;
        if (order == 1) return;
        double r_mm = ra2 * ir, r_aa = rm2 * ir, r_ma = -r_m * r_a * ir;
        double S_mm = r_mm * v0;
        double S_ma = r_ma * v0 + u_m * v1;
        double S_aa = r_aa * v0 + 2.0 * r_a * v1 + u * v2;
        double T_mm = b * S_mm, T_ma = b * S_ma, T_aa = b * S_aa;
        double P_mm = Q2 * T_m * T_m + Q1 * T_mm;
        double P_ma = Q2 * T_m * T_a + Q1 * T_ma;
        double P_aa = Q2 * T_a * T_a + Q1 * T_aa;
        w = only ? 0 : 3;
        v[w + 0] = -yy * im * im + P_mm;
        v[w + 1] = P_ma;
        v[w + 2] = yy * ia * ia + p[2] + P_aa;
        if (order == 2) return;
        double r_mmm = -3.0 * ra2 * r_m * ir2;
        double r_mma = (2.0 * r_a - 3.0 * ra2 * r_a) * ir2;
        double r_maa = (-r_m + 3.0 * r_m * ra2) * ir2;
        double r_aaa = -3.0 * rm2 * r_a * ir2;
        double S_mmm = r_mmm * v0;
        double S_mma = r_mma * v0 + r_mm * v1;
        double S_maa = r_maa * v0 + 2.0 * r_ma * v1 + u_m * v2;
        double S_aaa = r_aaa * v0 + 3.0 * r_aa * v1 + 3.0 * r_a * v2 +
            u * v3;
        double T_mmm = b * S_mmm, T_mma = b * S_mma, T_maa = b * S_maa,
            T_aaa = b * S_aaa;
        double P_mmm = Q3 * T_m * T_m * T_m + 3.0 * Q2 * T_mm * T_m +
            Q1 * T_mmm;
        double P_mma = Q3 * T_m * T_m * T_a +
            Q2 * (T_mm * T_a + 2.0 * T_ma * T_m) + Q1 * T_mma;
        double P_maa = Q3 * T_m * T_a * T_a +
            Q2 * (T_aa * T_m + 2.0 * T_ma * T_a) + Q1 * T_maa;
        double P_aaa = Q3 * T_a * T_a * T_a + 3.0 * Q2 * T_aa * T_a +
            Q1 * T_aaa;
        w = only ? 0 : 6;
        v[w + 0] = 2.0 * yy * im * im * im + P_mmm;
        v[w + 1] = P_mma;
        v[w + 2] = P_maa;
        v[w + 3] = -2.0 * yy * ia * ia * ia + p[3] + P_aaa;
        if (order == 3) return;
        double r_mmmm = (-3.0 * ra2 + 15.0 * ra2 * rm2) * ir3;
        double r_mmma = (-6.0 * r_a * r_m + 15.0 * ra2 * r_a * r_m) * ir3;
        double r_mmaa = (2.0 - 15.0 * ra2 + 15.0 * ra2 * ra2) * ir3;
        double r_maaa = (9.0 * r_m * r_a - 15.0 * r_m * ra2 * r_a) * ir3;
        double r_aaaa = (-3.0 * rm2 + 15.0 * rm2 * ra2) * ir3;
        double S_mmmm = r_mmmm * v0;
        double S_mmma = r_mmma * v0 + r_mmm * v1;
        double S_mmaa = r_mmaa * v0 + 2.0 * r_mma * v1 + r_mm * v2;
        double S_maaa = r_maaa * v0 + 3.0 * r_maa * v1 + 3.0 * r_ma * v2 +
            u_m * v3;
        double S_aaaa = r_aaaa * v0 + 4.0 * r_aaa * v1 + 6.0 * r_aa * v2 +
            4.0 * r_a * v3 + u * v4;
        double T_mmmm = b * S_mmmm, T_mmma = b * S_mmma,
            T_mmaa = b * S_mmaa, T_maaa = b * S_maaa, T_aaaa = b * S_aaaa;
        double P_mmmm = Q4 * T_m * T_m * T_m * T_m +
            6.0 * Q3 * T_mm * T_m * T_m +
            Q2 * (3.0 * T_mm * T_mm + 4.0 * T_mmm * T_m) + Q1 * T_mmmm;
        double P_mmma = Q4 * T_m * T_m * T_m * T_a +
            Q3 * (3.0 * T_ma * T_m * T_m + 3.0 * T_mm * T_m * T_a) +
            Q2 * (3.0 * T_mm * T_ma + 3.0 * T_mma * T_m + T_mmm * T_a) +
            Q1 * T_mmma;
        double P_mmaa = Q4 * T_m * T_m * T_a * T_a +
            Q3 * (T_mm * T_a * T_a + T_aa * T_m * T_m +
                  4.0 * T_ma * T_m * T_a) +
            Q2 * (T_mm * T_aa + 2.0 * T_ma * T_ma +
                  2.0 * T_maa * T_m + 2.0 * T_mma * T_a) +
            Q1 * T_mmaa;
        double P_maaa = Q4 * T_m * T_a * T_a * T_a +
            Q3 * (3.0 * T_ma * T_a * T_a + 3.0 * T_aa * T_m * T_a) +
            Q2 * (3.0 * T_aa * T_ma + 3.0 * T_maa * T_a + T_aaa * T_m) +
            Q1 * T_maaa;
        double P_aaaa = Q4 * T_a * T_a * T_a * T_a +
            6.0 * Q3 * T_aa * T_a * T_a +
            Q2 * (3.0 * T_aa * T_aa + 4.0 * T_aaa * T_a) + Q1 * T_aaaa;
        w = only ? 0 : 10;
        v[w + 0] = -6.0 * yy * im * im * im * im + P_mmmm;
        v[w + 1] = P_mmma;
        v[w + 2] = P_mmaa;
        v[w + 3] = P_maaa;
        v[w + 4] = 6.0 * yy * ia * ia * ia * ia + p[4] + P_aaaa;
}

inline void pig1_d5_row(double yy, double mu, double sg, double* v) {
    if (ISNAN(yy)) {
        for (int j = 0; j < 6; ++j) v[j] = yy;
        return;
    }
    if (!R_finite(yy) || yy < 0.0 || yy != std::floor(yy)) {
        for (int j = 0; j < 6; ++j) v[j] = R_NaN;
        return;
    }
    const double w = d7::pig1_w(mu, sg);
    double L[6];
    d7::pig5_logS_w(yy, w, L);
    const double L1 = L[1], L2 = L[2], L3 = L[3], L4 = L[4], L5 = L[5];
    double out_mu_mu_mu_mu_mu = 0.0, out_mu_mu_mu_mu_sigma = 0.0,
        out_mu_mu_mu_sigma_sigma = 0.0, out_mu_mu_sigma_sigma_sigma = 0.0,
        out_mu_sigma_sigma_sigma_sigma = 0.0,
        out_sigma_sigma_sigma_sigma_sigma = 0.0;
    double t0 = std::pow(mu, 5);
    double t1 = std::pow(sg, 6);
    double t2 = mu*sg;
    double t3 = 2*t2 + 1;
    double t4 = std::pow(t3, -11.0/2.0);
    double t5 = t1*t4;
    double t6 = (945.0/2.0)*L1;
    double t7 = std::pow(sg, 7);
    double t8 = std::pow(t3, -6);
    double t9 = t7*t8;
    double t10 = std::pow(sg, 8);
    double t11 = L3/std::pow(t3, 13.0/2.0);
    double t12 = L4/std::pow(t3, 7);
    double t13 = std::pow(sg, 5);
    double t14 = std::pow(t3, -5);
    double t15 = t13*t14;
    double t16 = 384*yy;
    double t17 = std::pow(t3, -7.0/2.0);
    double t18 = std::pow(sg, 4);
    double t19 = t17*t18;
    double t20 = std::sqrt(t3);
    double t21 = t20 + 1;
    double t22 = std::pow(t21, -2);
    double t23 = 150*t22;
    double t24 = std::pow(t3, -3);
    double t25 = t18*t24;
    double t26 = std::pow(t21, -3);
    double t27 = 300*t26;
    double t28 = std::pow(t3, -5.0/2.0);
    double t29 = std::pow(t21, -4);
    double t30 = t28*t29;
    double t31 = 360*t30;
    double t32 = std::pow(t21, -5);
    double t33 = std::pow(t3, -2);
    double t34 = t32*t33;
    double t35 = 240*t18;
    double t36 = std::pow(t3, -9.0/2.0);
    double t37 = std::pow(t3, -4);
    double t38 = std::pow(t21, -6);
    double t39 = (1.0/2.0)*L4;
    double t40 = std::pow(sg, 3);
    double t41 = t37*t40;
    double t42 = 192*yy;
    double t43 = t4*t6;
    double t44 = L2*t8;
    double t45 = t14*t18;
    double t46 = 72*t22;
    double t47 = std::pow(sg, 2);
    double t48 = t28*t47;
    double t49 = 144*t47;
    double t50 = t26*t33;
    double t51 = std::pow(t3, -3.0/2.0);
    double t52 = t29*t51;
    double t53 = mu*t17*t40;
    double t54 = 240*t22;
    double t55 = t24*t40;
    double t56 = 480*t26;
    double t57 = mu*t40;
    double t58 = 576*t30;
    double t59 = 384*t34;
    double t60 = std::pow(mu, 2);
    double t61 = t18*t36;
    double t62 = 210*t22;
    double t63 = t18*t37;
    double t64 = 420*t26;
    double t65 = 540*t29;
    double t66 = 480*t32;
    double t67 = t28*t60;
    double t68 = mu*t51;
    double t69 = -1.0/2.0*sg*t68 + (1.0/2.0)/t20;
    double t70 = L2*t69;
    double t71 = t13*t36;
    double t72 = L3*t69;
    double t73 = L4*t69;
    double t74 = sg*t24;
    double t75 = 48*yy;
    double t76 = t17*t47;
    double t77 = 90*L1;
    double t78 = t40*t60;
    double t79 = t14*t16;
    double t80 = t60*t76;
    double t81 = t24*t47;
    double t82 = t47*t60;
    double t83 = std::pow(mu, 3);
    double t84 = 30*t17;
    double t85 = t40*t84;
    double t86 = (45.0/4.0)*t63;
    double t87 = (3.0/4.0)*t71;
    double t88 = std::pow(t69, 2);
    double t89 = (1.0/8.0)*t1*t36;
    double t90 = sg*t67;
    double t91 = -t68 + (3.0/2.0)*t90;
    double t92 = L2*t91;
    double t93 = (15.0/2.0)*t19;
    double t94 = L3*t91;
    double t95 = (9.0/4.0)*t13*t37;
    double t96 = 9*t28;
    double t97 = t47*t83;
    double t98 = 144*t60;
    double t99 = t17*t83;
    double t100 = sg*t99;
    double t101 = sg*t83;
    double t102 = t36*t62;
    double t103 = std::pow(mu, 4);
    double t104 = t103*t47;
    double t105 = t37*t64;
    double t106 = 240*t103;
    double t107 = sg*t96;
    double t108 = 3*t81;
    double t109 = L4*t88;
    double t110 = 45*mu*t76;
    double t111 = 12*mu*t41;
    double t112 = std::pow(t69, 3);
    double t113 = (1.0/4.0)*t25;
    double t114 = (105.0/2.0)*t36;
    double t115 = t114*t78;
    double t116 = 12*t45*t60;
    double t117 = (3.0/2.0)*mu;
    double t118 = t117*t63;
    double t119 = 3*mu;
    double t120 = -15.0/2.0*t100 + (9.0/2.0)*t28*t60;
    double t121 = L2*t120;
    double t122 = (3.0/2.0)*t28*t40;
    double t123 = (15.0/2.0)*t53;
    double t124 = -t123 + (9.0/2.0)*t28*t47;
    double t125 = -2*t68 + 3*t90;
    double t126 = t125*t69;
    double t127 = -t118 + t122*t69 + t55;
    double t128 = (9.0/2.0)*t28;
    double t129 = t128*t47;
    double t130 = t108 - t111 + t116 + t122*t91 - t123*t69 + t129*t69;
    double t131 = t17*t60;
    double t132 = sg*t103;
    double t133 = sg*t0;
    double t134 = t17*t65;
    double t135 = 240*t28*t38;
    double t136 = mu*t96;
    double t137 = L3*t88;
    double t138 = 3*t92;
    double t139 = sg*t51;
    double t140 = L4*t112;
    double t141 = 45*sg*t131;
    double t142 = t2*t28;
    double t143 = (1.0/2.0)*t47;
    double t144 = t143*t51;
    double t145 = t114*t97;
    double t146 = t117*t48;
    double t147 = t114*t132 - t83*t84;
    double t148 = (15.0/2.0)*t80;
    double t149 = -t139 + t146;
    double t150 = 2*t139;
    double t151 = L3*t126;
    double t152 = -15*t100 + 9*t28*t60;
    double t153 = t51*t69;
    double t154 = t143*t153;
    double t155 = sg*t128*t60 - 3*t68;
    double t156 = t155*t88;
    double t157 = t125*t91;
    double t158 = 3*t70;
    double t159 = 6*mu*sg*t28 - t148 - t51;
    double t160 = -t139*t69 - t144*t91 + t146*t69 + t149*t69;
    double t161 = -6*mu*sg*t28*t69 - 3*mu*t28*t47*t91 + t120*t144 + t148*t69 - t149*t91 + t150*t91 + t153;
    double t162 = std::pow(mu, 6);
    double t163 = t126 + t69*t91;
    double t164 = t120*t69;
    double t165 = t152*t69 + t157 + std::pow(t91, 2);
    // mu_mu_mu_mu_mu
    out_mu_mu_mu_mu_mu = -975.0/4.0*L2*t9 - 1.0/32.0*L5*std::pow(sg, 10)/std::pow(t3, 15.0/2.0) + 540*mu*t13*t17*t29 + 210*mu*t13*t22*t36 + 480*mu*t13*t24*t32 + 420*mu*t13*t26*t37 + 240*mu*t13*t28*t38 - 15.0/8.0*std::pow(sg, 9)*t12 - 285.0/8.0*t10*t11 - t15*t16 - t18*t31 - t19*t23 - t25*t27 - t34*t35 - t5*t6 + 24*yy/t0;
    // mu_mu_mu_mu_sigma
    out_mu_mu_mu_mu_sigma = (525.0/2.0)*L1*t18*t36 + (261.0/2.0)*L2*t15 + (63.0/4.0)*L3*t5 + (1.0/16.0)*L5*t10*t69*t8 - 435.0/2.0*mu*t1*t44 - 3.0/4.0*mu*t10*t12 - 99.0/4.0*mu*t11*t7 - mu*t13*t43 - mu*t16*t45 - mu*t55*t56 + (87.0/4.0)*t1*t14*t72 + t19*t60*t65 + t25*t60*t66 + t35*t38*t67 + t39*t9 + (9.0/4.0)*t4*t7*t73 + t41*t42 + t46*t48 + t49*t50 + t49*t52 - t53*t54 - t57*t58 - t57*t59 + t60*t61*t62 + t60*t63*t64 + (105.0/2.0)*t70*t71;
    // mu_mu_mu_sigma_sigma
    out_mu_mu_mu_sigma_sigma = 420*L1*mu*t36*t40 + 180*L2*mu*t14*t18 + (105.0/2.0)*L2*mu*t18*t36*t69 - 45*L2*t41 + L2*t69*((105.0/2.0)*mu*t61 - t85) + 18*L3*mu*t13*t14*t69 + (27.0/2.0)*L3*mu*t13*t4 - 15.0/4.0*L3*t61 + L3*t69*(18*mu*t13*t14 - t69*t93 - t86) + (9.0/8.0)*L4*mu*t1*t4*t69 + L4*t69*((9.0/8.0)*mu*t1*t4 - t69*t95 - t87) - L4*t89*t91 - L5*t88*t89 + 108*mu*sg*t22*t28 + 216*mu*sg*t26*t33 + 216*mu*sg*t29*t51 + 288*mu*t37*t47*yy - 99.0/8.0*t1*t11*t60 - 180*t13*t44*t60 + 540*t17*t29*t40*t83 - t18*t43*t60 + 210*t22*t36*t40*t83 - 12*t22*t51 - 270*t22*t80 + 480*t24*t32*t40*t83 + 420*t26*t37*t40*t83 - 540*t26*t60*t81 - 24*t26/t3 + 240*t28*t38*t40*t83 - 648*t29*t47*t67 - 432*t34*t82 - t70*t85 - t72*t86 - t73*t87 - t74*t75 - t76*t77 - t78*t79 - t92*t93 - t94*t95;
    // mu_mu_sigma_sigma_sigma
    out_mu_mu_sigma_sigma_sigma = -135*L1*t17*t2 + L1*t96 - 54*L2*mu*t37*t47 + 144*L2*t14*t78 + 6*L2*t74 + L3*t113*t120 + L4*t113*t126 + L5*t112*t113 - mu*t24*t75 + 288*sg*t37*t60*yy - t100*t54 - t101*t58 - t101*t59 + t102*t104 + t103*t65*t76 + t103*t66*t81 + t104*t105 + t106*t38*t48 + t107*t70 + t108*t72 - t109*t118 + t109*t55 - t110*t70 - t111*t72 + t115*t70 + t116*t72 - t119*t63*t94 + t121*t122 + t124*t92 + t127*t94 + t130*t72 - 120*t18*t44*t83 + t36*t6*t82 - t40*t43*t83 + t46*t67 + t47*t92*t96 + t50*t98 + t52*t98 - 15*t53*t92 + 2*t55*t94 - t56*t74*t83 + 2*t70*(t107 - t110 + t115) + t72*(t124*t69 + t130) + t73*(t113*t91 - t118*t69 + t127*t69 + t55*t69) - t79*t97;
    // mu_sigma_sigma_sigma_sigma
    out_mu_sigma_sigma_sigma_sigma = 420*L1*t101*t36 - L2*t144*t147 - L3*t144*t157 - L3*t152*t154 - L5*t144*std::pow(t69, 4) + mu*t121*t129 + t0*t66*t74 + t102*t133 - t103*t17*t23 - t103*t24*t27 - t103*t31 - t104*t43 + t105*t133 - t106*t34 + t119*t151*t48 - 3*t121*t139 + t121*t149 - t131*t77 - t132*t79 + t133*t134 + t133*t135 + t136*t70 + 6*t137*t142 - t137*t148 - t137*t51 + t138*t159 - t138*t51 - t139*t140 + t140*t146 - t141*t70 + 18*t142*t92 + t145*t70 - t150*t151 - t156*t39*t47*t51 + t158*(t136 - t141 + t145) + t160*t94 + t37*t42*t83 + t72*(t159*t69 - t161) + t72*(2*t159*t69 - t161) + t73*(-t125*t154 - t139*t88 + t146*t88 + t160*t69) - 45.0/2.0*t80*t92;
    // sigma_sigma_sigma_sigma_sigma
    out_sigma_sigma_sigma_sigma_sigma = L1*((525.0/2.0)*t103*t36 - 945.0/2.0*t133*t4) + L5*std::pow(t69, 5) - t0*t79 + t102*t162 + t105*t162 + 2*t120*t92 + t121*t125 + t134*t162 + t135*t162 + t137*(-45.0/2.0*t100 + (27.0/2.0)*t28*t60) + t138*t152 + t140*(-4*t68 + 6*t90) + t147*t158 + t151*t155 + t162*t24*t66 + t163*t94 + t70*(105*t132*t36 - 60*t99) + t72*(t164 + t165) + t72*(2*t164 + t165) + t73*(t156 + t163*t69);
    v[0] = out_mu_mu_mu_mu_mu;
    v[1] = out_mu_mu_mu_mu_sigma;
    v[2] = out_mu_mu_mu_sigma_sigma;
    v[3] = out_mu_mu_sigma_sigma_sigma;
    v[4] = out_mu_sigma_sigma_sigma_sigma;
    v[5] = out_sigma_sigma_sigma_sigma_sigma;
}

inline void pig2_d5_row(double yy, double mu, double al, double* v) {
    if (ISNAN(yy)) {
        for (int j = 0; j < 6; ++j) v[j] = yy;
        return;
    }
    if (!R_finite(yy) || yy < 0.0 || yy != std::floor(yy)) {
        for (int j = 0; j < 6; ++j) v[j] = R_NaN;
        return;
    }
    const double w = 1.0 / (2.0 * al);
    double L[6];
    d7::pig5_logS_w(yy, w, L);
    const double L1 = L[1], L2 = L[2], L3 = L[3], L4 = L[4], L5 = L[5];
    const double rr = std::hypot(mu, al);
    double out_mu_mu_mu_mu_mu = 0.0, out_mu_mu_mu_mu_alpha = 0.0,
        out_mu_mu_mu_alpha_alpha = 0.0, out_mu_mu_alpha_alpha_alpha = 0.0,
        out_mu_alpha_alpha_alpha_alpha = 0.0,
        out_alpha_alpha_alpha_alpha_alpha = 0.0;
    double t0 = 1.0/rr;
    double t1 = std::pow(rr, -2);
    double t2 = mu + rr;
    double t3 = 1.0/t2;
    double t4 = al + rr;
    double t5 = std::pow(t4, -4);
    double t6 = t3*t5;
    double t7 = std::pow(mu, 2);
    double t8 = 6*t7;
    double t9 = al*t6*t8;
    double t10 = std::pow(t2, -2);
    double t11 = std::pow(t4, -3);
    double t12 = t10*t11;
    double t13 = std::pow(t2, -3);
    double t14 = std::pow(t4, -2);
    double t15 = t13*t14;
    double t16 = 2*t13;
    double t17 = t16*yy;
    double t18 = 1.0/t4;
    double t19 = mu*t18;
    double t20 = t19 + 1;
    double t21 = std::pow(t2, -4);
    double t22 = 6*t21;
    double t23 = t20*t22;
    double t24 = al*t23;
    double t25 = mu*t24 - t17;
    double t26 = al*t12*t8 + al*t15*t8 + t25 + t9;
    double t27 = std::pow(rr, -3);
    double t28 = t14*t3;
    double t29 = t10*t20;
    double t30 = al*t29;
    double t31 = al*t28*t7 + mu*t30 - t3*yy;
    double t32 = t11*t3;
    double t33 = 2*t32*t7;
    double t34 = al*t33;
    double t35 = 2*t10;
    double t36 = t35*t7;
    double t37 = t14*t36;
    double t38 = t16*t20;
    double t39 = al*t38;
    double t40 = mu*t39;
    double t41 = -al*t37 + t10*yy - t34 - t40;
    double t42 = 2*t1;
    double t43 = t41*t42;
    double t44 = 4*t32;
    double t45 = al*t44;
    double t46 = mu*t45;
    double t47 = t10*t14;
    double t48 = 4*t47;
    double t49 = al*t48;
    double t50 = mu*t49;
    double t51 = t16*t19;
    double t52 = -6*al*mu*t20*t21 - 2*al*t10*t11*t7 - 4*al*t13*t14*t7 + al*t51 + t17 + t39 + t46 + t50;
    double t53 = mu*t0*t26 + 2*mu*t27*t31 - mu*t43 - t52;
    double t54 = t0*t53;
    double t55 = t0*t41;
    double t56 = t1*t31;
    double t57 = t10*t19;
    double t58 = 2*t28;
    double t59 = al*t58;
    double t60 = t47*t7;
    double t61 = al*t60;
    double t62 = al*t57 + mu*t59 + t10*yy + t30 - t40 - t61;
    double t63 = mu*t55 - mu*t56 + t62;
    double t64 = t1*t63;
    double t65 = -t52;
    double t66 = mu*t1;
    double t67 = t14*t7;
    double t68 = al*t16*t67;
    double t69 = 4*t13;
    double t70 = t20*t69;
    double t71 = t19*t69;
    double t72 = al*t18*t35 - al*t70 - al*t71 + t25 - t50 + t59 + t68;
    double t73 = t55 - t56;
    double t74 = mu*t0*t65 - t62*t66 + t72 + t73;
    double t75 = mu*t54 - mu*t64 + t74;
    double t76 = t1*t75;
    double t77 = 24*al;
    double t78 = std::pow(t4, -5);
    double t79 = t3*t78;
    double t80 = t7*t79;
    double t81 = t77*t80;
    double t82 = t10*t5;
    double t83 = t7*t82;
    double t84 = t11*t7;
    double t85 = t13*t84;
    double t86 = t21*t67;
    double t87 = std::pow(t2, -5);
    double t88 = t20*t87;
    double t89 = t77*t88;
    double t90 = mu*t89;
    double t91 = -6*t21*yy + t90;
    double t92 = -t77*t83 - t77*t85 - t77*t86 - t81 - t91;
    double t93 = t0*t92;
    double t94 = std::pow(rr, -4);
    double t95 = 6*t94;
    double t96 = t31*t95;
    double t97 = 3*t66;
    double t98 = 6*t27;
    double t99 = t41*t98;
    double t100 = 18*t86;
    double t101 = 12*al;
    double t102 = t101*t85;
    double t103 = t19*t22;
    double t104 = t101*t12;
    double t105 = t101*t15;
    double t106 = t22*yy - t90;
    double t107 = t8*t82;
    double t108 = al*t107;
    double t109 = t101*t6;
    double t110 = mu*t109 - t108;
    double t111 = -al*t100 + al*t103 + mu*t104 + mu*t105 - t102 + t106 + t110 + t24;
    double t112 = mu*t93 - mu*t96 + mu*t99 + t111 - t26*t97;
    double t113 = t0*t112;
    double t114 = t42*t53;
    double t115 = 2*t27;
    double t116 = t115*t63;
    double t117 = mu*t0;
    double t118 = t18*t69;
    double t119 = al*t118;
    double t120 = mu*t42;
    double t121 = mu*t115;
    double t122 = t20*t21;
    double t123 = 12*t122;
    double t124 = t101*t86;
    double t125 = 8*t12;
    double t126 = al*t125;
    double t127 = mu*t126;
    double t128 = t19*t21;
    double t129 = 16*mu;
    double t130 = t129*t15;
    double t131 = t0*t26;
    double t132 = t115*t31;
    double t133 = t69*t84;
    double t134 = al*t133;
    double t135 = t131 + t132 - t134 - t43;
    double t136 = al*t123 + al*t130 + t101*t128 + t106 + t111*t117 - t119 - t120*t65 + t121*t62 - t124 + t127 + t135 - t45 - t49;
    double t137 = mu*t113 - mu*t114 + mu*t116 + t136;
    double t138 = -8*al*mu*t10*t11;
    double t139 = t134 + t138 + t45;
    double t140 = 16*al*mu*t13*t14 + 12*al*mu*t18*t21 + 12*al*t20*t21 - t119 - t124 - t139 - t49 - t91;
    double t141 = 6*t47;
    double t142 = t13*t18;
    double t143 = 12*t142;
    double t144 = t22*t67;
    double t145 = al*t144;
    double t146 = -12*al*mu*t13*t14 - 18*al*mu*t18*t21 + al*t141 + al*t143 - 18*al*t20*t21 + t145 + t91;
    double t147 = -mu*t0*t140 - 2*t0*t65 + t146 + t42*t62 + t66*t72;
    double t148 = -t0*t53 + t64;
    double t149 = -mu*t0*t136 + t147 + t148 + t66*t74;
    double t150 = mu*t0*t137 - mu*t76 - t149;
    double t151 = t0*t150;
    double t152 = -t149;
    double t153 = 2*t0;
    double t154 = -t147;
    double t155 = 3*t0;
    double t156 = -t146;
    double t157 = 24*yy;
    double t158 = 96*al;
    double t159 = 72*al;
    double t160 = t18*t21;
    double t161 = t19*t87;
    double t162 = t15*t77;
    double t163 = t157*t87;
    double t164 = std::pow(t2, -6);
    double t165 = al*mu;
    double t166 = -t163 + 120*t164*t165*t20;
    double t167 = t162 + t166;
    double t168 = al*t14;
    double t169 = mu*t21;
    double t170 = t168*t169;
    double t171 = 48*t170;
    double t172 = t67*t87;
    double t173 = t172*t77;
    double t174 = -t171 + t173;
    double t175 = -t158*t161 - t158*t88 + t159*t160 + t167 + t174;
    double t176 = 3*t1;
    double t177 = al*t160;
    double t178 = 48*al;
    double t179 = t14*t169;
    double t180 = t104 - t159*t179;
    double t181 = t21*t84;
    double t182 = t101*t181;
    double t183 = t11*t13;
    double t184 = mu*t77;
    double t185 = -t183*t184;
    double t186 = t182 + t185;
    double t187 = -t159*t161 - t159*t88 + t167 + t172*t178 + 36*t177 + t180 + t186;
    double t188 = t117*t187 + t140*t155 - t156*t66 + t175 - t176*t72;
    double t189 = 4*t1;
    double t190 = 4*t27;
    double t191 = al*t88;
    double t192 = 48*t165;
    double t193 = t183*t192;
    double t194 = al*t181;
    double t195 = 36*t194;
    double t196 = t13*t5*t7;
    double t197 = t101*t196;
    double t198 = -t184*t82;
    double t199 = t109 + t198;
    double t200 = t197 + t199;
    double t201 = t101*t160 + t105 + t159*t172 - t161*t178 + t166 + t180 - 48*t191 - t193 + t195 + t200;
    double t202 = t111*t153 + t117*t201 - t120*t140 + t121*t72 + t187 - t189*t65 + t190*t62;
    double t203 = t117*t202 + t136*t153 - t154*t66 + t188 - t42*t74;
    double t204 = t192*t79;
    double t205 = -120*al*mu*t164*t20 - 24*al*t10*t7*t78 - 72*al*t11*t21*t7 - 48*al*t13*t5*t7 - 96*al*t14*t7*t87 + t161*t77 + t163 + t171 + t192*t82 + t193 + t204 + t89;
    double t206 = mu*t95;
    double t207 = mu*t98;
    double t208 = t176*t26;
    double t209 = -t208 + t93 - t96 + t99;
    double t210 = -t111*t97 - t117*t205 + t201 - t206*t62 + t207*t65 + t209;
    double t211 = t113 - t114 + t116;
    double t212 = t117*t210 - t120*t136 + t121*t74 + t202 + t211;
    double t213 = t0*t137 - t76;
    double t214 = t117*t212 - t152*t66 + t203 + t213;
    double t215 = 120*al;
    double t216 = t215*t3*t7/std::pow(t4, 6);
    double t217 = t10*t7*t78;
    double t218 = t166 + t172*t215 + 120*t194 + t196*t215 + t215*t217 + t216;
    double t219 = std::pow(rr, -5);
    double t220 = 24*mu;
    double t221 = t41*t94;
    double t222 = mu*t0*t218 + 24*mu*t219*t31 + 12*mu*t26*t27 - t205 - t220*t221 - 4*t66*t92;
    double t223 = -t112*t97 + t117*t222 - t206*t63 + t207*t53 + t210;
    double t224 = t117*t223 - t120*t137 + t121*t75 + t212;
    double t225 = t117*t224 - t150*t66 + t214;
    double t226 = -t220*t88;
    double t227 = 4*t12;
    double t228 = -t37 + t9;
    double t229 = -mu*t38 + t68;
    double t230 = al*t227*t7 + t228 + t229 - t33;
    double t231 = t0*t230;
    double t232 = -mu*t10*t20 - t14*t3*t7 + t34 + t61;
    double t233 = -t232;
    double t234 = t1*t233;
    double t235 = t11*t36;
    double t236 = al*t235;
    double t237 = t168*t35;
    double t238 = -mu*t237 + mu*t58 + t229 + t236 + t29 - t46 + t57 - t60;
    double t239 = mu*t231 - mu*t234 + t238;
    double t240 = t1*t239;
    double t241 = 18*t83;
    double t242 = -6*t10*t11*t7 + t81;
    double t243 = -6*mu*t20*t21 + t145;
    double t244 = al*t241 + t102 - 6*t13*t14*t7 + t242 + t243 - 6*t3*t5*t7;
    double t245 = -t244;
    double t246 = t230*t42;
    double t247 = mu*t246;
    double t248 = -4*t13*t14*t7;
    double t249 = mu*t44;
    double t250 = 8*t85;
    double t251 = al*t250;
    double t252 = mu*t48;
    double t253 = -12*al*mu*t3*t5 + t108 + t252;
    double t254 = -4*al*mu*t13*t14 - 2*t10*t11*t7 + t138 + t243 + t248 + t249 + t251 + t253 + t38 + t51;
    double t255 = mu*t0*t245 + 2*mu*t233*t27 - t247 - t254;
    double t256 = -t254;
    double t257 = -8*al*mu*t13*t14 - 2*t10*t18 - 2*t13*t14*t7 + t139 - 2*t14*t3 + t237 + t243 + t252 + t70 + t71;
    double t258 = -mu*t0*t256 - t0*t230 + t234 + t238*t66 + t257;
    double t259 = -mu*t0*t255 + mu*t240 + t258;
    double t260 = -t259;
    double t261 = t0*t260;
    double t262 = -t258;
    double t263 = -t257;
    double t264 = al*t0;
    double t265 = 12*mu;
    double t266 = t173 + t226;
    double t267 = t15*t265 + t266;
    double t268 = t104 + t105 + 18*t122 + 18*t128 - t141 - t143 - t144 - 36*t170 + t186 + t267;
    double t269 = -12*t86;
    double t270 = t168*t69;
    double t271 = t130 + t200 - t48;
    double t272 = mu*t125 - t118 + t123 + t126 + 12*t128 - t133 - 32*t165*t183 - t179*t77 + t181*t77 + t266 + t269 + t270 + t271 - t44;
    double t273 = t117*t272 + t153*t256 - t238*t42 - t263*t66 + t268;
    double t274 = t0*t255;
    double t275 = t0*t245;
    double t276 = t115*t233;
    double t277 = -12*t85;
    double t278 = 36*al*t196;
    double t279 = t217*t77;
    double t280 = -t204 + t279;
    double t281 = t12*t265 + t280;
    double t282 = -t100 - t101*t179 + t103 - t107 - 36*t165*t82 + t185 + t195 + t23 + t265*t6 + t267 + t277 + t278 + t281;
    double t283 = t117*t282 - t120*t256 + t121*t238 - t246 + t272 + t275 + t276;
    double t284 = t117*t283 - t240 - t262*t66 + t273 + t274;
    double t285 = -24*t85;
    double t286 = 24*t7;
    double t287 = -t286*t82;
    double t288 = t216 + t287;
    double t289 = t158*t217 + t159*t196 + t178*t181 + t266 + t285 - t286*t79 + t288 - 24*t86;
    double t290 = t117*t289 - t206*t233 + t207*t230 - t245*t97 + t282;
    double t291 = t117*t290 - t120*t255 + t121*t239 + t283;
    double t292 = t117*t291 - t260*t66 + t284;
    double t293 = t228 + t236 - t44*t7;
    double t294 = t1*t293;
    double t295 = t101*t83;
    double t296 = t134 + t248;
    double t297 = 8*t10*t11*t7 - t295 - t296 + 12*t3*t5*t7 - t81;
    double t298 = 8*t32;
    double t299 = -4*al*mu*t10*t11 + mu*t298 - 4*t10*t11*t7 + t253 + t296;
    double t300 = -mu*t0*t297 + mu*t294 + t299;
    double t301 = -t300;
    double t302 = t0*t301;
    double t303 = -t299;
    double t304 = al*t1;
    double t305 = al*t183;
    double t306 = t182 + t269;
    double t307 = t12*t129 + t306;
    double t308 = al*t227 - t129*t305 - t250 + t271 - t298 + t307;
    double t309 = t0*t297;
    double t310 = 12*t7;
    double t311 = 8*mu;
    double t312 = t15*t311 + t196*t77 + t198 + t220*t6 + t280 - t305*t311 + t307 - t310*t82 - 16*t85;
    double t313 = t117*t312 - t294 - t303*t66 + t308 + t309;
    double t314 = t216 - 36*t83;
    double t315 = t159*t217 + t278 + t285 + t306 + t314 - 48*t80;
    double t316 = t117*t315 - t120*t297 + t121*t293 + t312;
    double t317 = t117*t316 - t301*t66 + t313;
    double t318 = t108 + t242 - 18*t3*t5*t7;
    double t319 = -t318;
    double t320 = t0*t319;
    double t321 = al*t0*t137 - al*t76 - t259;
    double t322 = t213 - t260*t304 + t264*t291 + t317;
    double t323 = al*t42;
    double t324 = al*t115;
    double t325 = mu*t6;
    double t326 = t197 + t277;
    double t327 = -mu*t101*t82 - t241 + t281 + 36*t325 + t326;
    double t328 = t178*t217 + t314 + t326 - 72*t80;
    double t329 = t117*t328 - t319*t66 + t327;
    double t330 = al*t54 - al*t64 + t239;
    double t331 = t1*t330;
    double t332 = t67*t69;
    double t333 = al*t113 - al*t114 + al*t116 + mu*t23 + mu*t270 + mu*t275 + mu*t276 + t110 + t127 - t145 + t235 - t247 - t249 - t251 - t252 + t332 - t38 - t51;
    double t334 = -al*t0*t255 + al*t240 + t148 + t300;
    double t335 = al*t0*t333 - al*t331 - t334;
    double t336 = -t334;
    double t337 = -2*t240 + t264*t316 + 2*t274 - t301*t304 + t329;
    double t338 = t211 + t239*t324 - t255*t323 + t264*t290 + t316;
    double t339 = t0*t333 + t264*t338 - t304*t336 - t331 + t337;
    double t340 = al*t95;
    double t341 = al*t176;
    double t342 = al*t98;
    double t343 = t279 + t288 - 96*t80;
    double t344 = t264*t328 - 3*t294 - t304*t319 + 3*t309 + t343;
    double t345 = -al*t0*t297 + al*t294 - 2*t0*t230 + 2*t234 + t318;
    double t346 = -t345;
    double t347 = al*t0*t41 - al*t56 - t232;
    double t348 = t1*t347;
    double t349 = al*t131 + al*t132 - al*t43 + t230;
    double t350 = al*t231 - al*t234 + t293 + t73;
    double t351 = -al*t246 + al*t275 + al*t276 + t125*t7 + t135 - t295 + t310*t6 + t332 - t81;
    double t352 = -al*t0*t351 - t0*t349 + t304*t350 + t345 + t348;
    double t353 = -t352;
    double t354 = -t189*t230 + t190*t233 + t264*t315 + 2*t275 + t293*t324 - t297*t323 + t328;
    double t355 = t153*t351 + t264*t354 - t304*t346 + t344 - t350*t42;
    double t356 = al*t0*t349 - al*t348 + t350;
    double t357 = t1*t356;
    double t358 = al*t0*t92 - al*t208 + 6*al*t27*t41 - al*t96 - t244;
    double t359 = t0*t358;
    double t360 = t349*t42;
    double t361 = t115*t347;
    double t362 = al*t359 - al*t360 + al*t361 + t351;
    double t363 = al*t0*t362 - al*t357 - t352;
    double t364 = t209 + t230*t342 - t233*t340 - t245*t341 + t264*t289 + t315;
    double t365 = t264*t364 - t323*t351 + t324*t350 + t354 + t359 - t360 + t361;
    double t366 = t0*t362 + t264*t365 - t304*t353 + t355 - t357;
    out_mu_mu_mu_mu_mu = 4*t0*t156 + t117*t175 + t117*t188 + t117*t203 + t117*t214 + t117*t225 + t151 + t152*t153 + t154*t155 - t161*t215 + t166 + 120*t177 - 120*t191 + t157/std::pow(mu, 5);
    out_mu_mu_mu_mu_alpha = t117*t268 + t117*t273 + t117*t284 + t117*t292 + 24*t122 + 24*t128 - 24*t142 + t153*t262 + t155*t263 + t162 + t174 + t225*t264 + t226 + t261;
    out_mu_mu_mu_alpha_alpha = t104 + t117*t308 + t117*t313 + t117*t317 + t15*t220 + t151 + t153*t303 + t185 + t264*t292 + t264*(-t150*t304 + t224*t264 + t292) + t302 + t306 - 12*t47;
    out_mu_mu_alpha_alpha_alpha = t0*t321 + t117*t327 + t117*t329 + t12*t220 + t199 + 2*t261 + t264*t317 + t264*t322 + t264*(t264*(-t137*t323 + t223*t264 + t291 + t324*t75) - t304*t321 + t322) - 12*t32 + t320 + t326;
    out_mu_alpha_alpha_alpha_alpha = t0*t335 + t117*t343 + t153*t336 + t264*t329 + t264*t337 + t264*t339 + t264*(t264*(t264*(-t112*t341 + t222*t264 + t290 - t340*t63 + t342*t53) - t323*t333 + t324*t330 + t338) - t304*t335 + t339) + t280 + t287 + 3*t302 + 48*t325;
    out_alpha_alpha_alpha_alpha_alpha = -60*L1/std::pow(al, 6) - 60*L2/std::pow(al, 7) - 15*L3/std::pow(al, 8) - 5.0/4.0*L4/std::pow(al, 9) - 1.0/32.0*L5/std::pow(al, 10) + t0*t363 + t153*t353 + t155*t346 + t216 + t264*t343 + t264*t344 + t264*t355 + t264*t366 + t264*(t264*(t264*(t264*(-al*t189*t92 + t101*t26*t27 + t218*t264 + t219*t31*t77 - t221*t77 + t289) - t340*t347 - t341*t358 + t342*t349 + t364) - t323*t362 + t324*t356 + t365) - t304*t363 + t366) + 4*t320 - 120*t80 + t157/std::pow(al, 5);
    v[0] = out_mu_mu_mu_mu_mu;
    v[1] = out_mu_mu_mu_mu_alpha;
    v[2] = out_mu_mu_mu_alpha_alpha;
    v[3] = out_mu_mu_alpha_alpha_alpha;
    v[4] = out_mu_alpha_alpha_alpha_alpha;
    v[5] = out_alpha_alpha_alpha_alpha_alpha;
}

}  // namespace d7pigrow

#endif
