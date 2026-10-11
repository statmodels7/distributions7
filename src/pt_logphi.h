#ifndef D7_PT_LOGPHI_H
#define D7_PT_LOGPHI_H

#include <Rmath.h>
#include <cmath>
#include "pt_enet_logz.h"

// The derivatives g_k = d^k log Phi(t) / dt^k, k = 1..n (n <= 5).
//
// For t >= -1 they come from the inverse Mills ratio R = phi/Phi, formed on
// the log scale, and the recursion R' = -R(t + R):
//   g1 = R, g2 = -R(t + R), g3 = -g2(t + R) - R(1 + g2),
//   g4 = -g3(t + 2R) - 2 g2(1 + g2), g5 = -g4(t + 2R) - 3 g3(1 + 2 g2).
// Below t = -1 the recursion cancels: t + R is of size 1/|t| while t and R
// are of size |t|, and each further order compounds it. There
// log Phi(t) = -x^2/2 - log(2 pi)/2 + L(x), x = -t, L the logarithm of the
// Mills ratio, and L^(k)(x) = (-1)^k kappa_k with kappa_k the cumulants of
// the density proportional to exp(-x s - s^2/2) on s > 0 (pt_enet_logz.h):
//   g1 = x + kappa_1, g2 = -1 + kappa_2, g_k = kappa_k for k >= 3,
// with no difference of large numbers formed.

namespace d7 {

inline void logphi_derivs(double t, int n, double* g) {
    if (t >= -1.0) {
        const double R = std::exp(R::dnorm4(t, 0.0, 1.0, 1) -
                                  R::pnorm5(t, 0.0, 1.0, 1, 1));
        g[1] = R;
        if (n >= 2) g[2] = -R * (t + R);
        if (n >= 3) g[3] = -g[2] * (t + R) - R * (1.0 + g[2]);
        if (n >= 4) g[4] = -g[3] * (t + 2.0 * R) - 2.0 * g[2] * (1.0 + g[2]);
        if (n >= 5) g[5] = -g[4] * (t + 2.0 * R) - 3.0 * g[3] * (1.0 + 2.0 * g[2]);
        return;
    }
    const double x = -t;
    double m[12], kap[12];
    enet_s_moments(x, n, m);
    enet_cumulants(m, n, kap);
    g[1] = x + kap[1];
    if (n >= 2) g[2] = -1.0 + kap[2];
    for (int k = 3; k <= n; k++) g[k] = kap[k];
}

} // namespace d7

#endif
