#ifndef D7_PT_GAMMA1_H
#define D7_PT_GAMMA1_H

#include <Rcpp.h>
#include "psi_diff.h"
#include "pt_parts.h"

// Gamma in the mean and the dispersion: one function per component and
// order for the quantities the scalar registry reads, called by the vector
// kernels in gamma1.cpp and by the registry in d7_ccallable.cpp. With
// s = 1/phi and z = y/mu, every derivative in phi is the derivative in s
// carried across by s' = -s^2, s'' = 2 s^3 (gamma1.cpp has the full
// derivation). The functions are arithmetic only: the polygamma remainders
// they read are computed by the caller and passed in,
//
//   f1 = [log s - psi(s)] + [log z - (z - 1)]   (gamma1_f1 below)
//   f2 = 1/s - psi'(s),  f3 = -1/s^2 - psi''(s)  (psi1_rest, psi2_rest)

namespace d7 {

// the score's two cancelling pairs, each written out (psi_diff.h)
inline double gamma1_f1(double y, double m, double s) {
    double z = y / m, w = (y - m) / m;   // w = z - 1 without the rounding of z
    return d7::psi_log_rest(s) + d7::psi_Ew2(z, w);
}

inline double gamma1_score_mu(double y, double m, double s) {
    return s * ((y - m) / m) / m;
}

inline double gamma1_score_phi(double s, double f1) {
    return f1 * (-s * s);
}

inline double gamma1_hess_mu_mu(double y, double m, double s) {
    double z = y / m;
    double m2 = m * m;
    return s * (1.0 - 2.0 * z) / m2;
}

inline double gamma1_hess_phi_phi(double s, double f1, double f2) {
    double s1 = -s * s, s2 = 2.0 * s * s * s;
    return f2 * s1 * s1 + f1 * s2;
}

inline double gamma1_expected_mu_mu(double m, double s) {
    double m2 = m * m;
    return -s / m2;
}

inline double gamma1_expected_phi_phi(double s, double f2) {
    double s1 = -s * s;
    return f2 * s1 * s1;
}

inline double gamma1_dexpected_mu_mu_mu(double m, double s) {
    double im = 1.0 / m, im2 = im * im;
    return 2.0 * s * im2 * im;
}

// d_phi E_pp = -s^2 q'(s), q' = f3 s^4 + 4 f2 s^3
inline double gamma1_dexpected_phi_phi_phi(double s, double f2, double f3) {
    double s2 = s * s, s3 = s2 * s, s4 = s2 * s2;
    double q1 = f3 * s4 + 4.0 * f2 * s3;
    return -s2 * q1;
}

// The off-diagonal components that only the joint registry reads, written
// from the expressions of gamma1.cpp (gamma1_hessian_cpp and
// gamma1_dexpected1_cpp). E_mp is zero identically, so the routers write 0.0
// for it and no function is needed.

// d2 l / d mu d phi, with s1 = -s^2
inline double gamma1_hess_mu_phi(double y, double m, double s) {
    double s1 = -s * s;
    return (y - m) / m / m * s1;
}

// d E_mm / d phi = s^2 / mu^2
inline double gamma1_dexpected_mu_mu_phi(double m, double s) {
    double s2 = s * s;
    double im = 1.0 / m, im2 = im * im;
    return s2 * im2;
}

inline void gamma1_score_curv(int k, double y, const double* th,
                              double* out) {
    double m = th[0], s = 1.0 / th[1];
    if (k == 0) {
        out[0] = gamma1_score_mu(y, m, s);
        out[1] = gamma1_hess_mu_mu(y, m, s);
    } else {
        double f1 = gamma1_f1(y, m, s);
        double f2 = d7::psi1_rest(s);
        out[0] = gamma1_score_phi(s, f1);
        out[1] = gamma1_hess_phi_phi(s, f1, f2);
    }
}

inline void gamma1_info_dinfo(int k, double y, const double* th,
                              double* out) {
    double m = th[0], s = 1.0 / th[1];
    if (k == 0) {
        out[0] = gamma1_expected_mu_mu(m, s);
        out[1] = gamma1_dexpected_mu_mu_mu(m, s);
    } else {
        double f2 = d7::psi1_rest(s), f3 = d7::psi2_rest(s);
        out[0] = gamma1_expected_phi_phi(s, f2);
        out[1] = gamma1_dexpected_phi_phi_phi(s, f2, f3);
    }
}

// The routers of the joint registry, one per order: every component at one
// observation, on the parameter scale (mu, phi). The Hessian and the expected
// information fill a full 2 x 2 matrix and the derivative of the expected
// information a 2 x 2 x 2 array, column-major, out[a + 2 b + 4 c] being the
// derivative of E[l_ab] in parameter c. Every entry is written.
inline void gamma1_joint_grad(double y, const double* th, double* out) {
    double m = th[0], s = 1.0 / th[1];
    out[0] = gamma1_score_mu(y, m, s);
    out[1] = gamma1_score_phi(s, gamma1_f1(y, m, s));
}

inline void gamma1_joint_hess(double y, const double* th, double* out) {
    double m = th[0], s = 1.0 / th[1];
    out[0] = gamma1_hess_mu_mu(y, m, s);
    out[3] = gamma1_hess_phi_phi(s, gamma1_f1(y, m, s), d7::psi1_rest(s));
    out[1] = out[2] = gamma1_hess_mu_phi(y, m, s);
}

inline void gamma1_joint_info(double y, const double* th, double* out) {
    double m = th[0], s = 1.0 / th[1];
    out[0] = gamma1_expected_mu_mu(m, s);
    out[3] = gamma1_expected_phi_phi(s, d7::psi1_rest(s));
    out[1] = out[2] = 0.0;
}

inline void gamma1_joint_dinfo(double y, const double* th, double* out) {
    double m = th[0], s = 1.0 / th[1];
    for (int j = 0; j < 8; ++j) out[j] = 0.0;
    out[0] = gamma1_dexpected_mu_mu_mu(m, s);
    out[4] = gamma1_dexpected_mu_mu_phi(m, s);
    out[7] = gamma1_dexpected_phi_phi_phi(s, d7::psi1_rest(s),
                                          d7::psi2_rest(s));
}

// The quantities of the third and fourth derivatives, moved here from
// gamma1.cpp so that its vector kernels and the parts route below call the
// same functions: the derivatives of l in s (f1 to f4), in mu, and mixed.
struct Gamma1Fs {
    double f1, f2, f3, f4;      // derivatives in s
    double m1, m2, m3, m4;      // derivatives in mu
    double ms1, mms1, mmms1;    // mixed, one s
};

// the parts in mu and the mixed ones, arithmetic at one observation
inline void gamma1_fs_mu(double y, double m, double s, Gamma1Fs& o) {
    double z = y / m, w = (y - m) / m;   // w = z - 1 without the rounding of z
    double m2 = m * m, m3 = m2 * m, m4 = m2 * m2;
    o.m1 = s * w / m;
    o.m2 = s * (1.0 - 2.0 * z) / m2;
    o.m3 = s * (6.0 * z - 2.0) / m3;
    o.m4 = s * (6.0 - 24.0 * z) / m4;
    o.ms1 = w / m;
    o.mms1 = (1.0 - 2.0 * z) / m2;
    o.mmms1 = (6.0 * z - 2.0) / m3;
}

inline Gamma1Fs gamma1_fs(double y, double m, double phi) {
    double s = 1.0 / phi;
    Gamma1Fs o;
    // [log s - psi(s)] + [log z - (z-1)]: two cancelling pairs, each
    // written out.  See psi_diff.h.
    o.f1 = gamma1_f1(y, m, s);
    o.f2 = d7::psi1_rest(s);
    o.f3 = d7::psi2_rest(s);
    o.f4 = d7::psi3_rest(s);
    gamma1_fs_mu(y, m, s, o);
    return o;
}

// The expected versions of the same quantities, at E[z] = 1 and
// E[log z] = digamma(s) - log s.
inline Gamma1Fs gamma1_fs_expected(double m, double phi) {
    double s = 1.0 / phi;
    double m2 = m * m, m3 = m2 * m, m4 = m2 * m2;
    Gamma1Fs o;
    o.f1 = 0.0;
    o.f2 = d7::psi1_rest(s);
    o.f3 = d7::psi2_rest(s);
    o.f4 = d7::psi3_rest(s);
    o.m1 = 0.0;
    o.m2 = -s / m2;
    o.m3 = 4.0 * s / m3;
    o.m4 = -18.0 * s / m4;
    o.ms1 = 0.0;
    o.mms1 = -1.0 / m2;
    o.mmms1 = 4.0 / m3;
    return o;
}

// the unique components of the third derivative: mu mu mu, mu mu phi,
// mu phi phi, phi phi phi
inline void gamma1_fs_d3(const Gamma1Fs &o, double s, double *out) {
    double s1 = -s * s, s2 = 2.0 * s * s * s, s3 = -6.0 * s * s * s * s;
    out[0] = o.m3;                                   // mu mu mu
    out[1] = o.mms1 * s1;                            // mu mu phi
    out[2] = o.ms1 * s2;                             // mu phi phi
    out[3] = o.f3 * s1 * s1 * s1 + 3.0 * o.f2 * s1 * s2 + o.f1 * s3;
}

// the unique components of the fourth derivative, in the same order
inline void gamma1_fs_d4(const Gamma1Fs &o, double s, double *out) {
    double s1 = -s * s, s2 = 2.0 * s * s * s, s3 = -6.0 * s * s * s * s,
           s4 = 24.0 * s * s * s * s * s;
    out[0] = o.m4;                                   // mu mu mu mu
    out[1] = o.mmms1 * s1;                           // mu mu mu phi
    out[2] = o.mms1 * s2;                            // mu mu phi phi
    out[3] = o.ms1 * s3;                             // mu phi phi phi
    out[4] = o.f4 * s1 * s1 * s1 * s1
           + 6.0 * o.f3 * s1 * s1 * s2
           + 3.0 * o.f2 * s2 * s2
           + 4.0 * o.f2 * s1 * s3
           + o.f1 * s4;
}

// The parts route (pt_parts.h). The buffer: [0] level, [1] y, [2] mu,
// [3] s = 1/phi, [4] to [7] f1 to f4, each filled when the level reaches its
// order. The expected parts leave y as NaN and f1 as 0, its expectation.
inline void gamma1_parts(double y, const double* th, int level, double* b) {
    const double m = th[0], s = 1.0 / th[1];
    b[0] = level; b[1] = y; b[2] = m; b[3] = s;
    b[4] = level >= 1 ? gamma1_f1(y, m, s) : R_NaN;
    b[5] = level >= 2 ? d7::psi1_rest(s) : R_NaN;
    b[6] = level >= 3 ? d7::psi2_rest(s) : R_NaN;
    b[7] = level >= 4 ? d7::psi3_rest(s) : R_NaN;
}

inline void gamma1_parts_expected(const double* th, int level, double* b) {
    const double m = th[0], s = 1.0 / th[1];
    b[0] = level; b[1] = R_NaN; b[2] = m; b[3] = s;
    b[4] = 0.0;
    b[5] = level >= 2 ? d7::psi1_rest(s) : R_NaN;
    b[6] = level >= 3 ? d7::psi2_rest(s) : R_NaN;
    b[7] = level >= 4 ? d7::psi3_rest(s) : R_NaN;
}

inline void gamma1_parts_grad(const double* b, double* out) {
    out[0] = gamma1_score_mu(b[1], b[2], b[3]);
    out[1] = gamma1_score_phi(b[3], b[4]);
}

inline void gamma1_parts_hess(const double* b, double* out) {
    out[0] = gamma1_hess_mu_mu(b[1], b[2], b[3]);
    out[3] = gamma1_hess_phi_phi(b[3], b[4], b[5]);
    out[1] = out[2] = gamma1_hess_mu_phi(b[1], b[2], b[3]);
}

// from the expected parts
inline void gamma1_parts_info(const double* b, double* out) {
    out[0] = gamma1_expected_mu_mu(b[2], b[3]);
    out[3] = gamma1_expected_phi_phi(b[3], b[5]);
    out[1] = out[2] = 0.0;
}

inline void gamma1_parts_d3(const double* b, double* out) {
    Gamma1Fs o;
    o.f1 = b[4]; o.f2 = b[5]; o.f3 = b[6]; o.f4 = b[7];
    gamma1_fs_mu(b[1], b[2], b[3], o);
    double u[4];
    gamma1_fs_d3(o, b[3], u);
    sym_fill(2, 3, u, out);
}

inline void gamma1_parts_d4(const double* b, double* out) {
    Gamma1Fs o;
    o.f1 = b[4]; o.f2 = b[5]; o.f3 = b[6]; o.f4 = b[7];
    gamma1_fs_mu(b[1], b[2], b[3], o);
    double u[5];
    gamma1_fs_d4(o, b[3], u);
    sym_fill(2, 4, u, out);
}

} // namespace d7

#endif
