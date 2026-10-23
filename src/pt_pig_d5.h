// Fifth derivative in w of log S_y(w), the Bessel polynomial of pt_pig.h,
// written for the fifth-order kernels of the Poisson-inverse Gaussian
// (d5_pig1.cpp and d5_pig2.cpp). It extends pig_logS_w() by one order, with
// the same two routes and the same switch between them. The sum and the
// Debye expansion are the forms of pt_pig.h; the fifth order is written out
// term by term (sympy's cse output for the Debye chain, pig5_debye_F()).

#ifndef D7_PT_PIG_D5_H
#define D7_PT_PIG_D5_H

#include "pt_pig.h"

namespace d7 {

// The sum, for 1 <= y < 21, to the fifth order in w. The same two forms as
// pig_logS_w_sum(), one order further.
//
// X = y (y - 1) w <= PIG_W_SWITCH: the Taylor coefficients of S(w + t) by the
// Horner scheme with derivatives, c_r = S^(r)(w)/r!, and the cumulants of the
// ratios N_r = r! c_r / c_0 = S^(r)/S:
//   L1 = N1, L2 = N2 - N1^2, L3 = N3 - 3 N1 N2 + 2 N1^3, ...
//   L5 = N5 - 5 N1 N4 - 10 N2 N3 + 20 N1^2 N3 + 30 N1 N2^2
//        - 60 N1^3 N2 + 24 N1^5.
// As in pig_logS_w_sum(), the error grows like y^(r-1) as w tends to zero.
//
// X > PIG_W_SWITCH: the cumulants kappa_j of k under the weights t_k / S,
// about the mean, with kappa_4 = mu_4 - 3 mu_2^2 and kappa_5 = mu_5 -
// 10 mu_2 mu_3 for the central moments mu_j. Since w^r d^r/dw^r is the
// falling factorial of the operator w d/dw, the signed Stirling numbers of
// the first kind give
//   w^5 L5 = kappa_5 - 10 kappa_4 + 35 kappa_3 - 50 kappa_2 + 24 kappa_1.
inline void pig5_logS_w_sum(double y, double w, double L[6]) {
    const int n = (int) y;
    double a[20];
    a[0] = 1.0;
    for (int k = 0; k < n - 1; ++k)
        a[k + 1] = a[k] * (y + k) * (y - k - 1.0) / (k + 1.0);
    if (y * (y - 1.0) * w <= PIG_W_SWITCH) {
        double dd[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        for (int k = n - 1; k >= 1; --k) {
            for (int r = 5; r >= 1; --r) dd[r] = dd[r] * w + dd[r - 1];
            dd[0] = dd[0] * w + a[k];
        }
        const double Sr = w * dd[0];
        for (int r = 5; r >= 1; --r) dd[r] = dd[r] * w + dd[r - 1];
        dd[0] = dd[0] * w + 1.0;
        const double c0 = dd[0];
        const double n1 = dd[1] / c0;
        const double n2 = 2.0 * dd[2] / c0;
        const double n3 = 6.0 * dd[3] / c0;
        const double n4 = 24.0 * dd[4] / c0;
        const double n5 = 120.0 * dd[5] / c0;
        L[0] = std::log1p(Sr);
        L[1] = n1;
        L[2] = n2 - n1 * n1;
        L[3] = n3 - 3.0 * n1 * n2 + 2.0 * n1 * n1 * n1;
        L[4] = n4 - 4.0 * n1 * n3 - 3.0 * n2 * n2 + 12.0 * n1 * n1 * n2 -
            6.0 * n1 * n1 * n1 * n1;
        L[5] = n5 - 5.0 * n1 * n4 - 10.0 * n2 * n3 + 20.0 * n1 * n1 * n3 +
            30.0 * n1 * n2 * n2 - 60.0 * n1 * n1 * n1 * n2 +
            24.0 * n1 * n1 * n1 * n1 * n1;
        return;
    }
    double t[20] = {0.0};
    if (w <= 1.0) {
        double pw = 1.0;
        for (int k = 0; k < n; ++k) { t[k] = a[k] * pw; pw *= w; }
    } else {
        const double v = 1.0 / w;
        double u = 1.0;
        for (int k = n - 1; k >= 0; --k) { t[k] = a[k] * u; u *= v; }
    }
    double Sr = 0.0, s1 = 0.0;
    for (int k = n - 1; k >= 1; --k) {
        Sr += t[k];
        s1 += k * t[k];
    }
    const double S = t[0] + Sr;
    const double k1 = s1 / S;
    double c2 = 0.0, c3 = 0.0, c4 = 0.0, c5 = 0.0;
    for (int k = 0; k < n; ++k) {
        const double d = k - k1, d2 = d * d;
        c2 += t[k] * d2;
        c3 += t[k] * d2 * d;
        c4 += t[k] * d2 * d2;
        c5 += t[k] * d2 * d2 * d;
    }
    c2 /= S; c3 /= S; c4 /= S; c5 /= S;
    const double k2 = c2, k3 = c3, k4 = c4 - 3.0 * c2 * c2,
        k5 = c5 - 10.0 * c2 * c3;
    L[0] = (w <= 1.0) ? std::log1p(Sr) : (n - 1) * std::log(w) + std::log(S);
    const double iw = 1.0 / w;
    L[1] = k1 * iw;
    L[2] = (k2 - k1) * iw * iw;
    L[3] = (k3 - 3.0 * k2 + 2.0 * k1) * iw * iw * iw;
    L[4] = (k4 - 6.0 * k3 + 11.0 * k2 - 6.0 * k1) * iw * iw * iw * iw;
    L[5] = (k5 - 10.0 * k4 + 35.0 * k3 - 50.0 * k2 + 24.0 * k1) *
        iw * iw * iw * iw * iw;
}

// The fifth-order chain of the Debye branch, F_k = d^k/dz^k of
// nu f(z) + h(z) + log P(c(z)) for k = 1..5, from c = z/r, d = 1/r, r =
// sqrt(1 + z^2), with dc/dz = d^3, dd/dz = -c d^2, f' = d/(1 + d),
// h' = -c d/2, and lc[m] = l_{m+1} the cumulants of log P at c. Generated
// symbolically (pig1/gen_debye5.py in the section's scratch folder) and
// written out as cse temporaries; L_k = (2 nu)^k F_k.
inline void pig5_debye_F(double c, double d, double nu, const double* lc,
                         double F[6]) {
    const double l1 = lc[0], l2 = lc[1], l3 = lc[2], l4 = lc[3], l5 = lc[4];
    const double q0 = (1.0/2.0)*c;
    const double q1 = std::pow(d, 3);
    const double q2 = d + 1;
    const double q3 = 1.0/q2;
    const double q4 = nu*q3;
    const double q5 = std::pow(d, 4);
    const double q6 = std::pow(d, 6);
    const double q7 = std::pow(d, 2);
    const double q8 = std::pow(q2, -2);
    const double q9 = nu*q8;
    const double q10 = -d*q9 + 3*l1*q7 + nu*q3 - q0;
    const double q11 = q10*q7;
    const double q12 = std::pow(d, 9);
    const double q13 = std::pow(d, 7);
    const double q14 = 3*q13;
    const double q15 = l2*q14;
    const double q16 = q0*q7 - q11;
    const double q17 = 2*q1;
    const double q18 = std::pow(d, 5);
    const double q19 = 6*l1;
    const double q20 = 2*d;
    const double q21 = nu/std::pow(q2, 3);
    const double q22 = d*q19 + q20*q21 - 2*q9;
    const double q23 = q22*q7;
    const double q24 = c*d;
    const double q25 = 2*q24;
    const double q26 = -c*q23 + 6*l2*q18 - q10*q25 - q17;
    const double q27 = q26*q7;
    const double q28 = 9*std::pow(d, 10);
    const double q29 = l3*q28;
    const double q30 = std::pow(c, 2);
    const double q31 = -q14 + 12*q18*q30;
    const double q32 = q1*q31;
    const double q33 = c*d - q10*q20 - q23;
    const double q34 = q33*q7;
    const double q35 = -c*q34 - q15 + q18 - q27;
    const double q36 = std::pow(d, 8);
    const double q37 = 21*l2*q6;
    const double q38 = 6*q7;
    const double q39 = 4*q22;
    const double q40 = nu/std::pow(q2, 4);
    const double q41 = -6*d*q40 + q19 + 6*q21;
    const double q42 = q41*q7;
    const double q43 = 2*c;
    const double q44 = -c*q42 + 30*l2*q5 - q10*q43 - q24*q39 - q38;
    const double q45 = q44*q7;
    const double q46 = -c*q37 - c*q45 + 9*l3*q36 + q1*q33 + 3*q16*q7 - q25*q26;
    const double q47 = q46*q7;
    const double q48 = 60*q30*q5 - 21*q6;
    const double q49 = c*q7;
    const double q50 = l2*q1;
    const double q51 = 2*c + 2*d*nu*q8 - d*q39 - l1*q38 - 2*q4 - q42;
    const double q52 = -q20*q26 - q25*q33 - q37 - q45 - q49*q51 + 5*q5;
    const double q53 = 3*q7;
    const double F1 = -d*q0 + d*q4 + l1*q1;
    const double F2 = -c*q11 + l2*q6 - 1.0/2.0*q5;
    const double F3 = -c*q15 - c*q27 + l3*q12 + q1*q16;
    const double F4 = -c*q29 - c*q47 + std::pow(d, 12)*l4 + l2*q32 + q1*q35;
    const double F5 = -18*c*std::pow(d, 13)*l4 + std::pow(d, 15)*l5 + l3*q1*(-q28 + 63*q30*q36 + q32) + q1*(24*c*l2*q36 + q1*(-c*q17 - 2*q34) - q29 - q47 - q49*q52) - q49*(-90*c*l3*q12 + 12*std::pow(d, 11)*l4 + l2*q31*q53 + q1*q52 - q25*q46 + q35*q53 + q48*q50 - q49*(-126*c*l2*q18 + 6*d*q16 + 72*l3*q13 + q1*q51 - 4*q24*q44 - q26*q43 + 6*q33*q7 - q49*(-6*c*q22 - 12*d + 120*l2*q1 - 6*q24*q41 - q49*(24*d*nu/std::pow(q2, 5) - 24*q40)))) + q50*(24*c*q36 - q48*q49);
    F[1] = F1;
    F[2] = F2;
    F[3] = F3;
    F[4] = F4;
    F[5] = F5;
}

// The uniform expansion, for y >= 21, to the fifth order in w. The same
// quantities as pig_logS_w_debye() (nu, z, c, d, P and its Horner
// coefficients), with the cumulants of P in c carried to the fifth order:
// N_j = P^(j)(c)/P(c) = j! q_j / p0 and the cumulants of log P,
//   l1 = N1, l2 = N2 - N1^2, l3 = N3 - 3 N1 N2 + 2 N1^3, ...
// L_k = (2 nu)^k F_k, the chain of pig5_debye_F().
inline void pig5_logS_w_debye(double y, double w, double L[6]) {
    const double nu = y - 0.5;
    const double z = 2.0 * nu * w;
    const double r = std::hypot(1.0, z);
    const double d = 1.0 / r;
    const double c = (z > 1e150) ? 1.0 : z * d;
    int K = PIG_KD;
    for (int k = 0; k <= PIG_KD; ++k)
        if (nu > PIG_UK_NU[k]) { K = k; break; }
    double P[3 * PIG_KD + 1];
    for (int e = 0; e <= 3 * K; ++e) P[e] = 0.0;
    const double inu = 1.0 / nu;
    double f = 1.0;
    for (int k = 0; k <= K; ++k) {
        const double* uk = PIG_U + PIG_UOFF[k];
        for (int e = 0; e <= 3 * k; ++e) P[e] += f * uk[e];
        f *= -inu;
    }
    // Taylor coefficients of P at c, q[j] = P^(j)(c)/j!, by Horner
    double q[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int e = 3 * K; e >= 0; --e) {
        for (int i = 5; i >= 1; --i) q[i] = q[i] * c + q[i - 1];
        q[0] = q[0] * c + P[e];
    }
    const double p0 = q[0];
    // the same Horner without the constant term, for log1p(p0 - 1)
    double acc = 0.0;
    for (int e = 3 * K; e >= 1; --e) acc = acc * c + P[e];
    const double q0 = acc * c;
    const double n1 = q[1] / p0;
    const double n2 = 2.0 * q[2] / p0;
    const double n3 = 6.0 * q[3] / p0;
    const double n4 = 24.0 * q[4] / p0;
    const double n5 = 120.0 * q[5] / p0;
    const double l1 = n1;
    const double l2 = n2 - n1 * n1;
    const double l3 = n3 - 3.0 * n1 * n2 + 2.0 * n1 * n1 * n1;
    const double l4 = n4 - 4.0 * n1 * n3 - 3.0 * n2 * n2 +
        12.0 * n1 * n1 * n2 - 6.0 * n1 * n1 * n1 * n1;
    const double l5 = n5 - 5.0 * n1 * n4 - 10.0 * n2 * n3 +
        20.0 * n1 * n1 * n3 + 30.0 * n1 * n2 * n2 - 60.0 * n1 * n1 * n1 * n2 +
        24.0 * n1 * n1 * n1 * n1 * n1;
    // asinh z and log(1 + z^2) without overflow of the square
    double asz, lq;
    if (z > 1e100) {
        const double lz = std::log(2.0 * nu) + std::log(w);
        asz = 0.6931471805599453 + lz;
        lq = 2.0 * lz;
    } else {
        asz = std::asinh(z);
        lq = std::log1p(z * z);
    }
    L[0] = nu * (asz - c / (1.0 + d)) - 0.25 * lq + std::log1p(q0);
    const double lc[5] = {l1, l2, l3, l4, l5};
    double F[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    pig5_debye_F(c, d, nu, lc, F);
    const double two = 2.0 * nu;
    double pw = 1.0;
    for (int k = 1; k <= 5; ++k) {
        pw *= two;
        L[k] = pw * F[k];
    }
}

// The Bessel log-polynomial and its first five derivatives in w = 1/(2 alpha),
// L[r] = d^r log S_y / dw^r, L[0] = log S_y. Zero for y < 1, where S_y = 1.
inline void pig5_logS_w(double y, double w, double L[6]) {
    for (int r = 0; r < 6; ++r) L[r] = 0.0;
    if (y < 1.0) return;
    if (y < PIG_Y_DEBYE) pig5_logS_w_sum(y, w, L);
    else pig5_logS_w_debye(y, w, L);
}

}  // namespace d7

#endif  // D7_PT_PIG_D5_H
