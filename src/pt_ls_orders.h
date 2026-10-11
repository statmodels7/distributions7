#ifndef D7_PT_LS_ORDERS_H
#define D7_PT_LS_ORDERS_H

#include <cmath>

// The derivatives of a location-scale log-density l = -log(s) + g(z),
// z = (y - mu)/s, at orders one to five. Differentiating A(z)/s^m with
// d/dmu [A/s^m] = -A'/s^(m+1) and d/ds [A/s^m] = -(z A' + m A)/s^(m+1)
// gives, for the component with j derivatives in mu and k = n - j in s,
//
//   s^n d^n l = [k == n] (-1)^n (n - 1)! + sum_i C[n][k][i] z^(i - j) g_i(z),
//
// with integer coefficients C (generated with sympy; ls_tables.py in the
// audit scripts). The family supplies term(p, i) = z^p g_i(z), p <= i, in a
// form that stays finite where z^p and g_i do not; out[k] receives the
// component with k derivatives in s, in the order of deriv_names().

namespace d7 {

static const double LS_C[6][6][6] = {
    {{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}},
    {{0, -1, 0, 0, 0, 0}, {0, -1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}},
    {{0, 0, 1, 0, 0, 0}, {0, 1, 1, 0, 0, 0}, {0, 2, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}},
    {{0, 0, 0, -1, 0, 0}, {0, 0, -2, -1, 0, 0}, {0, -2, -4, -1, 0, 0}, {0, -6, -6, -1, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}},
    {{0, 0, 0, 0, 1, 0}, {0, 0, 0, 3, 1, 0}, {0, 0, 6, 6, 1, 0}, {0, 6, 18, 9, 1, 0}, {0, 24, 36, 12, 1, 0}, {0, 0, 0, 0, 0, 0}},
    {{0, 0, 0, 0, 0, -1}, {0, 0, 0, 0, -4, -1}, {0, 0, 0, -12, -8, -1}, {0, 0, -24, -36, -12, -1}, {0, -24, -96, -72, -16, -1}, {0, -120, -240, -120, -20, -1}}
};

// (-1)^n (n - 1)!, the n-th derivative of -log(s) times s^n
inline double ls_log_term(int n) {
    double f = 1.0;
    for (int i = 2; i < n; i++) f *= i;
    return (n % 2) ? -f : f;
}

template <class Term>
inline void ls_order(int n, double s, Term term, double* out) {
    const double inv = 1.0 / s;
    double invn = 1.0;
    for (int i = 0; i < n; i++) invn *= inv;
    for (int k = 0; k <= n; k++) {
        const int j = n - k;
        double acc = (k == n) ? ls_log_term(n) : 0.0;
        for (int i = (j > 1 ? j : 1); i <= n; i++) {
            const double c = LS_C[n][k][i];
            if (c != 0.0) acc += c * term(i - j, i);
        }
        out[k] = acc * invn;
    }
}

} // namespace d7

#endif
