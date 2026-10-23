#ifndef D7_PT_PIG_MASS_H
#define D7_PT_PIG_MASS_H

#include <cmath>
#include <algorithm>
#include "pt_pig_debye.h"

// The log mass of the Poisson-inverse Gaussian at large counts, in both
// parametrizations. The assembly of pt_pig.h (pig1_slot(0, 0, ...) and
// pig2_value()) adds terms of size y log y, y log mu, y asinh(mu/alpha) and
// lgamma(y + 1), whose sum is of size log y: at y = 1e12 it has lost all but
// a few digits and at 9e15 none. The functions below carry out that
// cancellation in the algebra instead.
//
// With nu = y - 1/2, r = sqrt(mu^2 + alpha^2), s = sqrt(nu^2 + alpha^2) and
// the Debye expansion of log S_y(alpha) (psi_derivs_debye()),
//
//   l = y log mu - y asinh(mu/alpha) + (r - mu - alpha) + log S_y(alpha)
//       - lgamma(y + 1),
//   -alpha + log S_y(alpha) = -nu eta - (1/4) log(1 + nu^2/alpha^2) + log P(c),
//   nu eta = s + nu log(alpha / (nu + s)),
//
// with c = nu/s and P the Debye series of K_nu. Writing r = mu (1 + 2 d_m) and
// s = nu (1 + 2 d_n) where d_m = alpha^2 / (2 mu (r + mu)) and
// d_n = alpha^2 / (2 nu (s + nu)), and lgamma(y + 1) by Stirling's series,
// every term of size y log y or y cancels and
//
//   l = [-y log1p(d_m) + (r - mu)] + [nu log1p(d_n) - (s - nu)]
//       + 1/2 + nu log1p(-1/(2y)) - log y - stirlerr(y)
//       + log(alpha) - (1/2) log s + log P(c) - (1/2) log(4 pi).
//
// The two brackets are of size alpha^2 / mu and differ by a quantity that is
// quadratic in k = nu - mu, so they are combined, with
//
//   D = (d_m - d_n) / (1 + d_n)
//     = alpha^2 k (nu + mu) / (mu (u nu + v mu) (v + nu)),  u = r, v = s,
//
// into
//
//   Phi = -nu (log1p(D) - D) + N - (1/2) log1p(d_m),
//   N = -k^2 alpha^2 (mu + nu) (mu v + u v + alpha^2)
//       / ((v + nu)(u + mu)(u + v) mu (u nu + v mu)),
//
// where N = -nu D + k (1 - (mu + nu)/(u + v)) after the factorization of the
// numerator (verified symbolically) and has no linear term in k left. For
// |D| >= 1/2 (mu much smaller than nu, or the reverse) the split would
// reintroduce the cancellation and
//
//   Phi = -nu (log1p(d_m) - log1p(d_n)) + k (alpha^2/(u + mu) + alpha^2/(v + nu))
//         / (u + v) - (1/2) log1p(d_m)
//
// is used. mu, nu and alpha are divided by the largest of the three first (N is
// homogeneous of degree one, D and the d's of degree zero), so no square of a
// large quantity is formed.
//
// Range of use. y an integer >= PIG_LOGMASS_BIG_Y, mu > 0, alpha > 0; measured
// against 80 digits for y from 21 to 9e15, mu/y from 1e-12 to 1e12 and sigma
// from 1e-100 to 1e100: the relative error is at most 8e-15. The functions
// neither allocate nor call the R API.

namespace d7 {

// The smallest count at which the form is used: it needs the Debye expansion,
// which psi_derivs() and pig_logS_w() also switch to at 21, and it is at least
// as accurate as the assembly of pt_pig.h from there on.
const double PIG_LOGMASS_BIG_Y = 21.0;

// The range of mu / y over which the form was measured.
const double PIG_LOGMASS_BIG_RATIO = 1e12;

inline bool pig_logmass_big_ok(double y) {
    return y >= PIG_LOGMASS_BIG_Y;
}

inline bool pig_logmass_big_ok(double y, double mu) {
    return y >= PIG_LOGMASS_BIG_Y && mu > 0.0 &&
        mu >= y / PIG_LOGMASS_BIG_RATIO && mu <= y * PIG_LOGMASS_BIG_RATIO;
}

// log P(c) = log1p(q0), q0 = P(c) - 1, for the Debye series of K_nu at
// c = nu / sqrt(nu^2 + alpha^2): the loop of psi_derivs_debye().
inline double pig_mass_debye_q0(double nu, double c) {
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
    double acc = 0.0;
    for (int e = 3 * K; e >= 1; --e) acc = acc * c + P[e];
    return acc * c;
}

// Stirling's correction lgamma(y + 1) - ((y + 1/2) log y - y + log(2 pi)/2),
// five terms, absolute error below 6e-18 at y = 21 and falling as y^-11.
inline double pig_mass_stirlerr(double y) {
    const double iy = 1.0 / y, i2 = iy * iy;
    return iy * (1.0 / 12.0 - i2 * (1.0 / 360.0 - i2 * (1.0 / 1260.0 -
        i2 * (1.0 / 1680.0 - i2 * (1.0 / 1188.0)))));
}

// log1p(x) - x. Below |x| = 1/2 it is the series of 2 atanh(t) - 2t in
// t = x / (2 + x) plus the closed term 2 t - x = -x^2 / (2 + x), summed to
// full relative precision; above it the difference is formed directly.
inline double pig_mass_log1pmx(double x) {
    if (std::fabs(x) > 0.5) return std::log1p(x) - x;
    const double t = x / (2.0 + x), t2 = t * t;
    double term = t * t2;
    double sm = term / 3.0;
    double j = 5.0;
    for (int it = 0; it < 40; ++it) {
        term *= t2;
        const double add = term / j;
        sm += add;
        if (std::fabs(add) <= 1e-18 * std::fabs(sm)) break;
        j += 2.0;
    }
    return -x * x / (2.0 + x) + 2.0 * sm;
}

// The log mass in (y, mu, alpha), the common body of both parametrizations.
inline double pig_mass_big_alpha(double y, double m, double a) {
    const double nu = y - 0.5;
    const double L = std::max(m, std::max(nu, a));
    const double il = 1.0 / L;
    const double ms = m * il, ns = nu * il, as = a * il;
    const double u = std::hypot(ms, as), v = std::hypot(ns, as);
    const double a2 = as * as;
    const double um = u + ms, vn = v + ns, uv = u + v;
    const double den = ms * (u * ns + v * ms);
    const double dm = 0.5 * (as / ms) * (as / um);
    const double dn = 0.5 * (as / ns) * (as / vn);
    const double k = ns - ms;
    const double D = a2 * k * (ns + ms) / (den * vn);
    double Phi;
    if (std::fabs(D) < 0.5) {
        const double N = -k * k * a2 * (ms + ns) * (ms * v + u * v + a2) /
            (vn * um * uv * den);
        Phi = -nu * pig_mass_log1pmx(D) + N * L - 0.5 * std::log1p(dm);
    } else {
        const double lD = std::log1p(dm) - std::log1p(dn);
        Phi = -nu * lD + k * L * (a2 / um + a2 / vn) / uv -
            0.5 * std::log1p(dm);
    }
    const double s = v * L;
    const double cst = 0.5 + nu * std::log1p(-0.5 / y) - std::log(y) -
        pig_mass_stirlerr(y);
    const double q0 = pig_mass_debye_q0(nu, nu / s);
    // (1/2) log(4 pi)
    return Phi + cst + std::log(a) - 0.5 * std::log(s) + std::log1p(q0) -
        1.2655121234846454;
}

// pig2 (mu, alpha)
inline double pig2_logmass_big(double y, double mu, double alpha) {
    return pig_mass_big_alpha(y, mu, alpha);
}

// pig1 (mu, sigma), alpha = sqrt(1 + 2 sigma mu) / sigma
inline double pig1_logmass_big(double y, double mu, double sigma) {
    const double a = std::sqrt(1.0 + 2.0 * sigma * mu) / sigma;
    return pig_mass_big_alpha(y, mu, a);
}

}  // namespace d7

#endif
