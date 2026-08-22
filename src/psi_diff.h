#ifndef D7_PSI_DIFF_H
#define D7_PSI_DIFF_H

#include <Rcpp.h>
#include <cmath>

// DIFFERENCES OF psi AND psi' AT ARGUMENTS AN INTEGER APART.
//
// Several families here carry a boundary limit at which a shape or a
// dispersion runs to infinity and the family tends to a simpler one -- the
// negative binomial to the Poisson, the beta-binomial to the binomial, the
// Student t to the gaussian.  Every derivative in that parameter vanishes
// there, so it is written as a difference of terms that agree to leading
// order, and the direct form loses every digit it has well before the values
// a fit reaches.  Measured:
//
//   negbin, dl/dtheta    wrong by 1.0e-03 at theta = 1e6, by 4.4 at 1e7, and
//                        it CHANGES SIGN at 1e8; a fit of 2000 counts drawn
//                        at a true theta of 100 reports 1.6e+07
//   betabinom, dl/dA     wrong by 5.6e-05 at S = 1e6, EXACTLY ZERO at 1e9
//                        where the value is 1.9e-17, and 1.9e+08 out at 1e12;
//                        a fit at a true concentration of 3000 reports
//                        1.7e+08, and on binomial data 3.1e+09
//
// The two functions below are what those rewrites share.  Each takes the
// SHIFT k, which is a count or a size and therefore an integer, and returns
// the difference with its own leading behaviour subtracted, so that a caller
// pairing it with that behaviour cancels symbolically rather than in double
// precision.  Above a measured crossover each is its asymptotic series, from
//
//   psi(x)  = log x - 1/(2x) - 1/(12x^2) + 1/(120x^4) - ...
//   psi'(x) = 1/x + 1/(2x^2) + 1/(6x^3) - 1/(30x^5) + ...
//
// with every difference of powers written so that no power of x is formed on
// its own -- and with the products a*b avoided, since they overflow past
// x = 1.3e154 while the values a log link can produce reach 1.8e308.

namespace d7 {

// psi(x + k) - psi(x) - log1p(k/x)
//   = k/(2ab) + k(2a+k)/(12 a^2 b^2) - k(2a+k)(a^2+b^2)/(120 a^4 b^4) + ...
// with a = x, b = x + k.  Zero at k = 0, by either branch.
inline double psi_A_rest(double k, double x) {
    if (x >= 100.0) {
        const double a = x, b = x + k;
        const double p = (k / a) / b;          // not k/(a*b): a*b overflows
        const double q = 2.0 / b + p;          // = (2a + k)/(a b)
        const double ia2 = 1.0 / (a * a), ib2 = 1.0 / (b * b);
        return p * (0.5 + q * (1.0 / 12.0 - (ia2 + ib2) / 120.0));
    }
    return R::digamma(x + k) - R::digamma(x) - std::log1p(k / x);
}

// psi'(x + k) - psi'(x) + k/(x(x + k))
//   = -k(2a+k)/(2 a^2 b^2) - k(a^2+ab+b^2)/(6 a^3 b^3)
//     + k(a^4+a^3b+a^2b^2+ab^3+b^4)/(30 a^5 b^5) - ...
inline double psi_T_rest(double k, double x) {
    if (x >= 100.0) {
        const double a = x, b = x + k;
        const double p = (k / a) / b;
        const double q = 2.0 / b + p;
        const double ia = 1.0 / a, ib = 1.0 / b;
        const double s2 = ia * ia + ia * ib + ib * ib;
        const double s4 = ia * ia * ia * ia + ia * ia * ia * ib +
            ia * ia * ib * ib + ia * ib * ib * ib + ib * ib * ib * ib;
        return -p * (q * 0.5 + s2 / 6.0 - s4 / 30.0);
    }
    return R::trigamma(x + k) - R::trigamma(x) + k / (x * (x + k));
}

// log1p(w) - w = -w^2/2 + w^3/3 - w^4/4 + w^5/5 - ...
inline double psi_Ew(double w) {
    if (std::fabs(w) < 1e-3) {
        return w * w * (-0.5 + w * (1.0 / 3.0 + w * (-0.25 + w * 0.2)));
    }
    return std::log1p(w) - w;
}


// A POLYGAMMA MINUS ITS OWN LEADING ASYMPTOTE, which is the shape the gamma
// carries four times over.  From
//
//   psi^(n)(z) ~ (-1)^(n-1) [ (n-1)!/z^n + n!/(2 z^(n+1))
//                             + sum_k B_2k (2k+n-1)!/((2k)! z^(2k+n)) ]
//
// the leading term cancels and what is left is one order down, so the direct
// difference loses digits as the shape grows: measured, log(s) - psi(s) is
// wrong by 1.3e-09 at s = 1e5, by 2.7e-04 at 1e11 and reads EXACTLY ZERO at
// 1e14 where the value is 5e-15, and the three derivative versions go the
// same way.
//
// ⚠️ Unlike the two families above, a gamma fit does NOT reach there in the
// ordinary course: a dispersion of 1e-4 gives s = 1.0e+04, where the loss is
// 1e-11.  Reaching 1e8 would take a coefficient of variation of 1e-4, which
// is a degenerate fit -- which is precisely when the derivatives should not
// be noise, and why these are written even though nothing routine needs them.
//
// The crossover is 50 for all four, measured: the series is within 1.2e-12
// there and the direct form still has all its digits.

// log(s) - psi(s) = 1/(2s) + 1/(12 s^2) - 1/(120 s^4) + 1/(252 s^6) - ...
inline double psi_log_rest(double s) {
    if (s >= 50.0) {
        const double u = 1.0 / s, u2 = u * u;
        return u * (0.5 + u * (1.0 / 12.0 + u2 * (-1.0 / 120.0 + u2 / 252.0)));
    }
    return std::log(s) - R::digamma(s);
}

// 1/s - psi'(s) = -1/(2 s^2) - 1/(6 s^3) + 1/(30 s^5) - 1/(42 s^7) + ...
inline double psi1_rest(double s) {
    if (s >= 50.0) {
        const double u = 1.0 / s, u2 = u * u;
        return u2 * (-0.5 + u * (-1.0 / 6.0 + u2 * (1.0 / 30.0 - u2 / 42.0)));
    }
    return 1.0 / s - R::trigamma(s);
}

// -1/s^2 - psi''(s) = 1/s^3 + 1/(2 s^4) - 1/(6 s^6) + 1/(6 s^8) - ...
inline double psi2_rest(double s) {
    if (s >= 50.0) {
        const double u = 1.0 / s, u2 = u * u;
        return u2 * u * (1.0 + u * (0.5 + u2 * (-1.0 / 6.0 + u2 / 6.0)));
    }
    return -1.0 / (s * s) - R::psigamma(s, 2);
}

// 2/s^3 - psi'''(s) = -3/s^4 - 2/s^5 + 1/s^7 - (4/3)/s^9 + ...
inline double psi3_rest(double s) {
    if (s >= 50.0) {
        const double u = 1.0 / s, u2 = u * u;
        return u2 * u2 * (-3.0 + u * (-2.0 + u2 * (1.0 - u2 * 4.0 / 3.0)));
    }
    return 2.0 / (s * s * s) - R::psigamma(s, 3);
}

}  // namespace d7

#endif  // D7_PSI_DIFF_H
