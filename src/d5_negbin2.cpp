#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "d7_par.h"
using namespace Rcpp;

// Fifth-order derivatives of the NB2 log-mass in (mu, theta), written out.
//
// With s = theta + mu the log-mass, up to a constant, is
//   l = lgamma(y + theta) - lgamma(theta) + theta log theta
//       - theta log s + y log mu - y log s.
// The lgamma part depends on theta alone, so it contributes only to the
// component of five theta derivatives, where it is psi^(4)(y+theta) -
// psi^(4)(theta). Every other component is rational in (mu, theta, y).
//
// The rational components are written as the factored forms that
// derive.py (sympy) produced, with the common denominators kept whole.
//
// The pure-theta component carries the cancellation that deriv4 carries one
// order lower. The elementary part is -24 y / s^5 plus a term in mu and
// theta alone, and the polygamma difference is +24 y / theta^5 to leading
// order, so the two cancel against each other as theta grows. Two forms are
// used, as in negbin_hd.cpp:
//   - theta > 100 max(y, 1) (or y < 1.5): psi^(4)(y+theta) - psi^(4)(theta)
//     is 24 sum_{j<y} 1/(theta+j)^5 exactly, and y = sum_{j<y} 1 lets the
//     -24 y / s^5 term join the sum. The summand is
//       24 (mu - j) Q_j / (t^5 s^5),  t = theta + j,
//       Q_j = s^4 + s^3 t + s^2 t^2 + s t^3 + t^4,
//     since 1/t^5 - 1/s^5 = (s - t) Q / (t^5 s^5) with s - t = mu - j, which
//     has no cancellation. The sum has y terms and is short in this branch.
//   - otherwise the difference is R::psigamma(y+theta, 4) - R::psigamma(theta, 4)
//     and the elementary part is written in full, as deriv4 does.
// The cut is the one deriv4 uses, th > 100 * y, with y < 1.5 always taking
// the sum (at most one term).

static inline void negbin2_d5_row(double y, double mu, double th, double* out) {
    const double s = th + mu;
    const double s2 = s * s, s3 = s2 * s, s4 = s2 * s2, s5 = s4 * s;
    const double mu2 = mu * mu, mu3 = mu2 * mu, mu4 = mu2 * mu2, mu5 = mu4 * mu;
    const double th2 = th * th, th3 = th2 * th, th4 = th2 * th2;

    // components with no derivative in theta beyond the rational part
    out[0] = 24.0 * th * (y * (5.0 * mu4 + 10.0 * mu3 * th + 10.0 * mu2 * th2
                               + 5.0 * mu * th3 + th4) - mu5) / (mu5 * s5);
    out[1] = -6.0 * (-mu + 3.0 * th + 4.0 * y) / s5;
    out[2] = -12.0 * (-mu + th + 2.0 * y) / s5;
    out[3] = -6.0 * (-3.0 * mu + th + 4.0 * y) / s5;
    out[4] = -24.0 * (y - mu) / s5;

    // the rest of the pure-theta component, which carries no y through a
    // polygamma: -6 (mu^5 + 5 mu^4 th + 10 mu^3 th^2 + 10 mu^2 th^3) / (th^4 s^5)
    const double rest = -6.0 * (mu5 + 5.0 * mu4 * th + 10.0 * mu3 * th2
                                + 10.0 * mu2 * th3) / (th4 * s5);

    if (y < 1.5 || th > 100.0 * y) {
        double acc = 0.0;
        for (double j = 0.0; j < y; j += 1.0) {
            const double t = th + j, t2 = t * t, t3 = t2 * t, t5 = t3 * t2;
            const double Q = s4 + s3 * t + s2 * t2 + s * t3 + t2 * t2;
            acc += (mu - j) * Q / (t5 * s5);
        }
        out[5] = 24.0 * acc + rest;
    } else {
        out[5] = R::psigamma(y + th, 4.0) - R::psigamma(th, 4.0)
            + rest - 24.0 * y / s5;
    }
}

// [[Rcpp::export]]
List negbin2_deriv5_cpp(NumericVector y, NumericVector mu, NumericVector theta,
                        int threads = 1) {
    const int n = std::max(y.size(), std::max(mu.size(), theta.size()));
    NumericVector m_mmmmm(n), m_mmmmt(n), m_mmmtt(n), m_mmttt(n), m_mtttt(n), m_tttttt(n);
    const int ny = y.size(), nm = mu.size(), nt = theta.size();

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t idx) {
        const int i = (int) idx;
        double out[6];
        negbin2_d5_row(y[i % ny], mu[i % nm], theta[i % nt], out);
        m_mmmmm[i] = out[0];
        m_mmmmt[i] = out[1];
        m_mmmtt[i] = out[2];
        m_mmttt[i] = out[3];
        m_mtttt[i] = out[4];
        m_tttttt[i] = out[5];
    });

    return List::create(
        Named("mu_mu_mu_mu_mu") = m_mmmmm,
        Named("mu_mu_mu_mu_theta") = m_mmmmt,
        Named("mu_mu_mu_theta_theta") = m_mmmtt,
        Named("mu_mu_theta_theta_theta") = m_mmttt,
        Named("mu_theta_theta_theta_theta") = m_mtttt,
        Named("theta_theta_theta_theta_theta") = m_tttttt
    );
}
