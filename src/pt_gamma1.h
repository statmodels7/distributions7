#ifndef D7_PT_GAMMA1_H
#define D7_PT_GAMMA1_H

#include <Rcpp.h>
#include "psi_diff.h"

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

} // namespace d7

#endif
