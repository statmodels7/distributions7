#ifndef D7_PT_WRAPPERS_H
#define D7_PT_WRAPPERS_H

#include <Rcpp.h>
#include <cmath>
#include <algorithm>

// The zero wrappers' components for the scalar registry, from the inner
// family's quantities, which the registry's router reads from the inner
// family's own entries and passes: f0 the parent's mass at zero (as the
// exponential of its log-mass), s0 and H0 the parent's (k, k) score and
// second derivative at zero, g and h the same at the observation, Ep and Dp
// the parent's (k, k) expected second derivative and its derivative in
// parameter k. Each expression is the one the R method writes
// (zero_inflated.R, zero_adjusted.R, zero_wrapper_dexpected.R), operation
// for operation, so that the two routes agree to the last bit. z is the
// wrapper's probability, zi or za.

namespace d7 {

// ---- zero_inflated() --------------------------------------------------------

inline double zi_l0(double z, double f0) { return z + (1 - z) * f0; }

inline double zi_score_parent(double y, double z, double f0, double g) {
    double w = (y == 0) ? ((1 - z) * f0) / zi_l0(z, f0) : 1;
    return w * g;
}
inline double zi_score_zi(double y, double z, double f0) {
    return (y == 0) ? (1 - f0) / zi_l0(z, f0) : -1 / (1 - z);
}
inline double zi_curv_parent(double y, double z, double f0, double s0,
                             double H0, double h) {
    if (y != 0) return h;
    double w0 = ((1 - z) * f0) / zi_l0(z, f0);
    return w0 * H0 + w0 * (1 - w0) * s0 * s0;
}
inline double zi_curv_zi(double y, double z, double f0) {
    if (y == 0) {
        double l0 = zi_l0(z, f0), a = 1 - f0;
        return -(a * a) / (l0 * l0);
    }
    double b = 1 - z;
    return -1 / (b * b);
}
inline double zi_info_parent(double z, double f0, double s0, double H0,
                             double Ep) {
    double l0 = zi_l0(z, f0);
    double w0 = ((1 - z) * f0) / l0;
    double h_zi_0 = w0 * H0 + w0 * (1 - w0) * s0 * s0;
    double contrib_pos = Ep - H0 * f0;
    return l0 * h_zi_0 + (1 - z) * contrib_pos;
}
inline double zi_info_zi(double z, double f0) {
    double l0 = zi_l0(z, f0);
    double s0z = (1 - f0) / l0, spz = -1 / (1 - z);
    return -(l0 * (s0z * s0z) + (1 - l0) * (spz * spz));
}
inline double zi_dinfo_parent(double z, double f0, double s0, double H0,
                              double Dp) {
    double l0 = zi_l0(z, f0);
    return (1 - z) * Dp +
        z * (1 - z) * f0 * (H0 * s0 + s0 * H0) / l0 +
        (z * z) * (1 - z) * f0 * s0 * s0 * s0 / (l0 * l0);
}
inline double zi_dinfo_zi(double z, double f0) {
    double l0 = zi_l0(z, f0), b = 1 - z;
    return R_pow(1 - f0, 3.0) / (l0 * l0) - (1 - f0) / (b * b);
}
inline double zi_logpdf(double y, double z, double lf) {
    double la = std::log1p(-z) + lf;
    if (y != 0) return la;
    double lz = std::log(z);
    double m = (lz >= la) ? lz : la;
    return m + std::log(std::exp(lz - m) + std::exp(la - m));
}

// The joint registry's zero-inflated components, for a pair (a, b) of the
// parent's parameters and the probability: s0a, s0b and H0ab the parent's
// score and Hessian at zero, Eab and Dabc the parent's expected Hessian and
// its derivative in parameter c. The expressions are those of the R methods
// (zero_inflated.R) and of zi_dexpected_entry() (zero_wrapper_dexpected.R);
// with a = b = c they are the (k, k) functions above.
inline double zi_hess_pair0(double z, double f0, double s0a, double s0b,
                            double H0ab) {
    double w0 = ((1 - z) * f0) / zi_l0(z, f0);
    return w0 * H0ab + w0 * (1 - w0) * s0a * s0b;
}
inline double zi_hess_cross0(double z, double f0, double s0a) {
    double l0 = zi_l0(z, f0);
    return (-f0 / (l0 * l0)) * s0a;
}
inline double zi_info_pair(double z, double f0, double s0a, double s0b,
                           double H0ab, double Eab) {
    double l0 = zi_l0(z, f0);
    double w0 = ((1 - z) * f0) / l0;
    double h_zi_0 = w0 * H0ab + w0 * (1 - w0) * s0a * s0b;
    double contrib_pos = Eab - H0ab * f0;
    return l0 * h_zi_0 + (1 - z) * contrib_pos;
}
inline double zi_info_cross(double z, double f0, double s0a) {
    return -(f0 / zi_l0(z, f0)) * s0a;
}
// d E_ab / d theta_c, a, b and c parameters of the parent
inline double zi_dinfo_ppp(double z, double f0, double s0a, double s0b,
                           double s0c, double H0ac, double H0bc, double Dabc) {
    double l0 = zi_l0(z, f0);
    return (1 - z) * Dabc +
        z * (1 - z) * f0 * (H0ac * s0b + s0a * H0bc) / l0 +
        (z * z) * (1 - z) * f0 * s0a * s0b * s0c / (l0 * l0);
}
// d E_ab / d z
inline double zi_dinfo_ppz(double z, double f0, double s0a, double s0b,
                           double Eab) {
    double l0 = zi_l0(z, f0);
    return -Eab + f0 * s0a * s0b *
        ((1 - 2 * z) / l0 - z * (1 - z) * (1 - f0) / (l0 * l0));
}
// d E_az / d theta_c
inline double zi_dinfo_pzp(double z, double f0, double s0a, double s0c,
                           double H0ac) {
    double l0 = zi_l0(z, f0);
    return -(f0 / l0) * (s0c * s0a + H0ac) +
        (1 - z) * (f0 * f0) * s0a * s0c / (l0 * l0);
}
// d E_az / d z
inline double zi_dinfo_pzz(double z, double f0, double s0a) {
    double l0 = zi_l0(z, f0);
    return f0 * s0a * (1 - f0) / (l0 * l0);
}
// d E_zz / d theta_c
inline double zi_dinfo_zzp(double z, double f0, double s0c) {
    double l0 = zi_l0(z, f0), a = 1 - f0;
    return f0 * s0c * (2 * a / l0 + (1 - z) * (a * a) / (l0 * l0) +
                       1 / (1 - z));
}

// ---- zero_adjusted(), discrete parent ---------------------------------------

inline double zad_score_parent(double y, double f0, double s0, double g) {
    if (y == 0) return 0;
    double correction = f0 / (1 - f0);
    return g + correction * s0;
}
inline double za_score_za(double y, double z) {
    return (y == 0) ? 1 / z : -1 / (1 - z);
}
inline double zad_hess_correction(double f0, double s0, double H0) {
    double denom = 1 - f0;
    double f_prime = f0 * s0;
    double f_second = f0 * (H0 + s0 * s0);
    return (denom * f_second + f_prime * f_prime) / (denom * denom);
}
inline double zad_curv_parent(double y, double f0, double s0, double H0,
                              double h) {
    if (y == 0) return 0;
    return h + zad_hess_correction(f0, s0, H0);
}
inline double za_curv_za(double y, double z) {
    if (y == 0) return -1 / (z * z);
    double b = 1 - z;
    return -1 / (b * b);
}
inline double zad_info_parent(double z, double f0, double s0, double H0,
                              double Ep) {
    double denom = 1 - f0;
    double E_trunc = (Ep - f0 * H0) / denom;
    return (1 - z) * (E_trunc + zad_hess_correction(f0, s0, H0));
}
inline double za_info_za(double z) { return -1 / (z * (1 - z)); }
inline double zad_dinfo_parent(double z, double f0, double s0, double H0,
                               double Ep, double Dp) {
    double q = 1 - f0;
    return (1 - z) * (Dp / q + Ep * f0 * s0 / (q * q) +
                      f0 * (s0 * s0 * s0 + H0 * s0 + s0 * H0) / (q * q) +
                      2 * (f0 * f0) * s0 * s0 * s0 / R_pow(q, 3.0));
}
inline double za_dinfo_za(double z) {
    double x = z * (1 - z);
    return (1 - 2 * z) / (x * x);
}
// The joint registry's zero-adjusted components, for a pair (a, b) of the
// parent's parameters, from the R methods (zero_adjusted.R) and
// za_disc_dexpected_entry() (zero_wrapper_dexpected.R); with a = b = c they
// are the (k, k) functions above. The blocks between the parent's parameters
// and the probability are zero.
inline double zad_hess_correction_pair(double f0, double s0a, double s0b,
                                       double H0ab) {
    double denom = 1 - f0;
    double f_second = f0 * (H0ab + s0a * s0b);
    return (denom * f_second + (f0 * s0a) * (f0 * s0b)) / (denom * denom);
}
inline double zad_info_pair(double z, double f0, double s0a, double s0b,
                            double H0ab, double Eab) {
    double denom = 1 - f0;
    double E_trunc = (Eab - f0 * H0ab) / denom;
    return (1 - z) * (E_trunc + zad_hess_correction_pair(f0, s0a, s0b, H0ab));
}
// d E_ab / d theta_c
inline double zad_dinfo_ppp(double z, double f0, double s0a, double s0b,
                            double s0c, double H0ac, double H0bc, double Eab,
                            double Dabc) {
    double q = 1 - f0;
    return (1 - z) * (Dabc / q + Eab * f0 * s0c / (q * q) +
                      f0 * (s0c * s0a * s0b + H0ac * s0b + s0a * H0bc) / (q * q) +
                      2 * (f0 * f0) * s0a * s0b * s0c / R_pow(q, 3.0));
}
// d E_ab / d z
inline double zad_dinfo_ppz(double f0, double s0a, double s0b, double Eab) {
    double q = 1 - f0;
    return -(Eab / q + f0 * s0a * s0b / (q * q));
}
inline double zad_logpdf(double y, double z, double f0, double lf) {
    if (y == 0) return std::log(z);
    return std::log(1 - z) + lf - std::log1p(-f0);
}

// ---- zero_adjusted(), continuous parent -------------------------------------

inline double zac_score_parent(double y, double g) { return (y == 0) ? 0 : g; }
inline double zac_curv_parent(double y, double h) { return (y == 0) ? 0 : h; }
inline double zac_info_parent(double z, double Ep) { return (1 - z) * Ep; }
inline double zac_dinfo_parent(double z, double Dp) { return (1 - z) * Dp; }
inline double zac_logpdf(double y, double z, double lf) {
    if (y == 0) return std::log(z);
    return std::log(1 - z) + lf;
}

// ---- folded() ----------------------------------------------------------------
//
// From the parent's log-density at y and -y (lp, lm) and its (k, k) score
// and second derivative there (gp, hp at y, gm, hm at -y): the weight
// w = f(y)/(f(y) + f(-y)) of fold_parts(), formed as it is there from the two
// logarithms (w = plogis(lp - lm), log L = max + log1p(exp(-|lp - lm|))), so
// that neither depends on the densities being representable; where both
// are zero the weight is one half and L is zero. The score w gp + (1 - w) gm and
// the second derivative R2 - R1^2 of log_deriv() over the block ratios
// R1 = w gp + (1 - w) gm, R2 = w (hp + gp^2) + (1 - w)(hm + gm^2).

struct FoldW { double fp, fm, L, logL, w; };

inline FoldW fold_w(double lp, double lm) {
    FoldW W;
    if (ISNAN(lp)) lp = R_NegInf;
    if (ISNAN(lm)) lm = R_NegInf;
    W.fp = std::exp(lp);
    W.fm = std::exp(lm);
    const double d = lp - lm;
    if (ISNAN(d)) {
        W.w = 0.5;
        W.logL = R_NegInf;
    } else {
        W.w = R::plogis(d, 0.0, 1.0, 1, 0);
        W.logL = std::max(lp, lm) + std::log1p(std::exp(-std::fabs(d)));
    }
    W.L = std::exp(W.logL);
    return W;
}
inline double fold_score(double w, double gp, double gm) {
    return w * gp + (1 - w) * gm;
}
inline double fold_curv(double w, double gp, double gm, double hp, double hm) {
    double r1 = w * gp + (1 - w) * gm;
    double r2 = w * (hp + gp * gp) + (1 - w) * (hm + gm * gm);
    return r2 - r1 * r1;
}
inline double fold_logpdf(double y, const FoldW& W) {
    if (y < 0) return R_NegInf;
    return W.logL;
}

} // namespace d7

#endif
