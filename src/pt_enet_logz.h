#ifndef D7_PT_ENET_LOGZ_H
#define D7_PT_ENET_LOGZ_H

#include <Rmath.h>
#include <cmath>

// The normalizing constant of the elastic net and its partial derivatives in
// the two rates a = lambda alpha and c = lambda (1 - alpha), to order four.
//
// With x = a/sqrt(c), Z = 2 M_0(x)/sqrt(c), where
//   M_k(x) = int_0^inf s^k exp(-x s - s^2/2) ds,
// so log M_0 = L, the logarithm of the Mills ratio. Writing the derivatives
// of L through G = x - 1/M_0 cancels for large x, both terms being of size x
// while G is of size 1/x, and every further derivative compounds it. Here
// they are written through the moments m_k = M_k/M_0 of the density
// proportional to exp(-x s - s^2/2) on s > 0, whose cumulants are
// L^(k)(x) = (-1)^k kappa_k: no difference of large numbers is formed.
//
// The moments come from the ratios r_k = M_k/M_(k-1), which satisfy
// r_k = k/(x + r_(k+1)) (integration by parts, M_(k+1) = k M_(k-1) - x M_k);
// the recursion is run downward from r_(K+1) = 0, where it is stable, for
// x >= 1, and the forward recursion from M_0 is used below 1.
//
// Two representations of log Z are used, each where it does not cancel:
//   x < 2:  log Z = log 2 + L(x) - log(c)/2, x = a c^(-1/2);
//   x >= 2: log Z = log 2 - log a + F(t),     t = c/a^2,
// where F(t) = log(x M_0(x)) and F^(j)(t) is the j-th cumulant of
// W = -U^2/2 under the density proportional to exp(-u - t u^2/2), u > 0,
// whose raw moments are E[W^j] = (-1/2)^j x^(2j) m_(2j). In both, the
// partial derivatives of a function of a monomial in (a, c) are taken with
// the Euler operator: for Phi = sum_m b_m v^m H^(j + m)(v), v the monomial,
// one differentiation in the other rate maps b to a combination of b, so the
// result is a short sum of the terms v^m H^(j + m)(v).

namespace d7 {

// m_0 .. m_n of the density exp(-x s - s^2/2) on s > 0, x > 0
inline void enet_s_moments(double x, int n, double* m) {
    m[0] = 1.0;
    if (x < 1.0) {
        // M_0 = Phi(-x)/phi(x), and M_(k+1) = k M_(k-1) - x M_k forward
        const double M0 = std::exp(R::pnorm(-x, 0.0, 1.0, 1, 1) -
                                   R::dnorm(x, 0.0, 1.0, 1));
        double Mm = M0, Mk = 1.0 - x * M0;
        if (n >= 1) m[1] = Mk / M0;
        for (int k = 1; k < n; k++) {
            const double Mn = k * Mm - x * Mk;
            Mm = Mk;
            Mk = Mn;
            m[k + 1] = Mk / M0;
        }
        return;
    }
    // the continued fraction converges slowly near x = 1
    const int K = x < 2.0 ? 2000 : (x < 10.0 ? 300 : 80);
    double r = 0.0, rr[12];
    for (int k = K; k >= 1; k--) {
        r = k / (x + r);
        if (k <= n) rr[k] = r;
    }
    for (int k = 1; k <= n; k++) m[k] = m[k - 1] * rr[k];
}

// cumulants kappa_1 .. kappa_n from raw moments mu_0 = 1 .. mu_n
inline void enet_cumulants(const double* mu, int n, double* kap) {
    for (int m = 1; m <= n; m++) {
        double s = mu[m];
        double binom = 1.0;  // choose(m - 1, i - 1)
        for (int i = 1; i < m; i++) {
            s -= binom * kap[i] * mu[m - i];
            binom = binom * (m - i) / i;
        }
        kap[m] = s;
    }
}

// L^(k)(x) = d^k log M_0 / dx^k, k = 1..n (n <= 5)
inline void enet_logM_derivs(double x, int n, double* L) {
    double m[12], kap[12];
    enet_s_moments(x, n, m);
    enet_cumulants(m, n, kap);
    for (int k = 1; k <= n; k++) L[k] = (k % 2) ? -kap[k] : kap[k];
}

// apply one Euler step to b (length len): b'_m = -h (j + m) b_m - q b_m,
// b'_(m + 1) += -h b_m
inline int enet_euler_step(double* b, int len, int j, double h, int q) {
    double nb[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    for (int m = 0; m < len; m++) {
        nb[m] += (-h * (j + m) - q) * b[m];
        nb[m + 1] += -h * b[m];
    }
    for (int m = 0; m <= len; m++) b[m] = nb[m];
    return len + 1;
}

// d[i][j] = d^(i + j) log Z / da^i dc^j for 1 <= i + j <= n, n <= 4
inline void enet_logz_ac(double a, double c, int n, double d[5][5]) {
    const double x = a / std::sqrt(c);
    if (x >= 2.0) {
        const double t = 1.0 / (x * x);
        double m[12], mu[6], F[6] = {0, 0, 0, 0, 0, 0};
        enet_s_moments(x, 2 * n, m);
        mu[0] = 1.0;
        double h = 1.0, x2 = 1.0;
        for (int j = 1; j <= n; j++) {
            h *= -0.5;
            x2 *= x * x;
            mu[j] = h * x2 * m[2 * j];
        }
        enet_cumulants(mu, n, F);
        // T(j, m) = t^m F^(j + m)(t)
        for (int ord = 1; ord <= n; ord++) {
            for (int j = 0; j <= ord; j++) {
                const int i = ord - j;
                double b[8] = {1, 0, 0, 0, 0, 0, 0, 0};
                int len = 1;
                for (int q = 0; q < i; q++) len = enet_euler_step(b, len, j, 2.0, q);
                double v = 0.0, tm = 1.0;
                for (int mm = 0; mm < len; mm++) {
                    if (b[mm] != 0.0) v += b[mm] * tm * F[j + mm];
                    tm *= t;
                }
                v *= std::pow(a, -2.0 * j - i);
                if (j == 0) {
                    double f = 1.0;
                    for (int k = 2; k < i; k++) f *= k;
                    v += ((i % 2) ? -f : f) * std::pow(a, -i);
                }
                d[i][j] = v;
            }
        }
        return;
    }
    double L[6];
    enet_logM_derivs(x, n, L);
    for (int ord = 1; ord <= n; ord++) {
        for (int j = 0; j <= ord; j++) {
            const int i = ord - j;
            double b[8] = {1, 0, 0, 0, 0, 0, 0, 0};
            int len = 1;
            for (int q = 0; q < j; q++) len = enet_euler_step(b, len, i, 0.5, q);
            // X(i, m) = x^m L^(i + m)(x); at i = 0 the m = 0 term is L itself,
            // whose coefficient is zero once j >= 1
            double v = 0.0, xm = 1.0;
            for (int mm = 0; mm < len; mm++) {
                if (i + mm >= 1 && b[mm] != 0.0) v += b[mm] * xm * L[i + mm];
                xm *= x;
            }
            v *= std::pow(c, -0.5 * i - j);
            if (i == 0) {
                double f = 1.0;
                for (int k = 2; k < j; k++) f *= k;
                // d^j/dc^j of -log(c)/2
                v += -0.5 * ((j % 2) ? f : -f) * std::pow(c, -j);
            }
            d[i][j] = v;
        }
    }
}

} // namespace d7

#endif
