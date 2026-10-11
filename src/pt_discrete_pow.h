#ifndef D7_PT_DISCRETE_POW_H
#define D7_PT_DISCRETE_POW_H

#include <cmath>

// Helpers for the derivatives in the mean of the count families, whose
// orders are sums of terms c / d^k.

namespace d7 {

// c / d^k by repeated division, so that d^k is never formed (it underflows
// for d below about 1e-308^(1/k) while c / d^k is still representable or
// overflows to its true infinity); a zero coefficient gives zero, which is
// the limit of the term and not the 0/0 that c / d^k would give.
inline double coef_over_pow(double c, double d, int k) {
    if (c == 0.0) return 0.0;
    double r = c;
    for (int j = 0; j < k; ++j) r /= d;
    return r;
}

// T_k = y / mu^k - (y + c) / (mu + c)^k, the bracket of the k-th derivative
// in the mean of the geometric (c = 1) and of the negative binomial with
// size c. It is written as y (1 - q^k) / mu^k - c / (mu + c)^k with
// q = mu / (mu + c) and log q = -log1p(c / mu). The two terms of the original
// form agree to a relative c / mu where y is near mu, which loses every digit
// for a large mean or a small size; here 1 - q^k = -expm1(k log q) carries
// that difference exactly and the remaining subtraction loses at most a
// factor k.
inline double count_bracket(double y, double m, double c, int k) {
    double one_mq = -std::expm1(-k * std::log1p(c / m));
    return coef_over_pow(y * one_mq, m, k) - coef_over_pow(c, m + c, k);
}

// U_j = mu^-j - (mu + c)^-j = (1 - q^j) / mu^j, the bracket of the expected
// orders, E[T_{j+1}] under E[y] = mu, with the same device.
inline double count_bracket_expected(double m, double c, int j) {
    double one_mq = -std::expm1(-j * std::log1p(c / m));
    return coef_over_pow(one_mq, m, j);
}

} // namespace d7

#endif
