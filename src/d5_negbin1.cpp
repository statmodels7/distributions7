#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cmath>
using namespace Rcpp;

// Fifth-order derivatives of the NB1 log-mass in (mu, theta), written out.
//
// The NB1 log-mass is l = G(r) + r B(theta) + C(theta) with r = mu/theta,
// G(r) = lgamma(y+r) - lgamma(r), B = -log(1+theta), C = y log theta -
// y log(1+theta). Two routes are used, as in negbin1_components() of
// R/negbin1_higher.R, and each observation takes its own:
//
//   theta >= 1: the size route. Every derivative of G is a polygamma
//     difference, G^(m)(r) = psi^(m-1)(y+r) - psi^(m-1)(r), written here as
//     D_{m-1}. The difference is taken by d5nb1_psi_shift(), which is
//     psi_shift_diff() of R/psi_diff.R: direct below r = 100 and the
//     asymptotic series above it, so that the two polygammas, which agree
//     to leading order for a large r, are not subtracted in double precision.
//     The coefficients of D_k are the ones derive.py (sympy) produced; the
//     elementary parts carry Lg = log1p(theta) and are the factored forms.
//
//   theta < 1: the form in mu + i theta. The sum of lgamma(y+r) - lgamma(r)
//     is sum_{i<y} log(r + i), and the y log(theta) it carries cancels
//     exactly against C(theta), which leaves
//       l = sum_{i<y} log(mu + i theta) - (mu/theta) log1p(theta)
//           - y log1p(theta) + const.
//     Every component of order five is then 24 S_b, with
//       S_b = sum_{i<y} i^b / (mu + i theta)^5,
//     corrected by -mu M_5 - 24 y/(1+theta)^5 in the pure-theta component
//     and by -M_4 in the component with one mu, where M_b is the b-th
//     derivative of log1p(theta)/theta. M_b uses the power series below
//     theta = 0.5 (80 terms, the one deriv4 uses) and the recursion
//       M_{m+1} = [ (-1)^m m!/(1+theta)^(m+1) - (m+1) M_m ] / theta
//     above it. The sum needs y terms, which deriv4 also pays.
//
// The cut is theta < 1, per observation, as nb1_exact_cut() in R.
//
// The six components follow deriv_names(c("mu", "theta"), 5): component k
// has 5 - k derivatives in mu and k in theta, so k = 0 is mu_mu_mu_mu_mu
// and k = 5 is theta_theta_theta_theta_theta.

// psi_shift_diff(n, k, x) of R/psi_diff.R: psi^(n)(x + k) - psi^(n)(x)
static inline double d5nb1_psi_shift(int n, double k, double x) {
    if (x < 100.0) {
        if (n == 0) return R::digamma(x + k) - R::digamma(x);
        return R::psigamma(x + k, (double) n) - R::psigamma(x, (double) n);
    }
    const double a = x;
    const double l1 = std::log1p(k / a);
    // 1/(a+k)^p - 1/a^p, with no power of a formed on its own where it would overflow
    auto dp = [&](double p) { return std::expm1(-p * l1) / std::pow(a, p); };
    if (n == 0) {
        return l1 - 0.5 * dp(1.0) - dp(2.0) / 12.0 + dp(4.0) / 120.0 - dp(6.0) / 252.0;
    }
    const double nd = (double) n;
    const double sgn = (n - 1) % 2 == 0 ? 1.0 : -1.0;              // (-1)^(n-1)
    const double f1 = std::tgamma(nd), f2 = std::tgamma(nd + 1.0);  // (n-1)!, n!
    const double f3 = std::tgamma(nd + 2.0);                        // (n+1)!
    const double f5 = std::tgamma(nd + 4.0);                        // (n+3)!
    const double f7 = std::tgamma(nd + 6.0);                        // (n+5)!
    return sgn * (f1 * dp(nd) + f2 / 2.0 * dp(nd + 1.0) + f3 / 12.0 * dp(nd + 2.0)
                  - f5 / 720.0 * dp(nd + 4.0) + f7 / 30240.0 * dp(nd + 6.0));
}

// M_b = d^b/dtheta^b [log1p(theta)/theta], b = 0..5
static inline void d5nb1_M(double th, double* M) {
    if (th < 0.5) {
        for (int b = 0; b <= 5; ++b) {
            double acc = 0.0, p = 1.0;
            for (int m = 0; m <= 80; ++m) {
                double co = 1.0;
                for (int i = 1; i <= b; ++i) co *= (m + i);
                const double sgn = (m % 2 == 0) ? 1.0 : -1.0;
                acc += sgn * co * p / (m + b + 1.0);
                p *= th;
            }
            M[b] = ((b % 2 == 0) ? 1.0 : -1.0) * acc;
        }
    } else {
        const double om = 1.0 + th;
        M[0] = std::log1p(th) / th;
        double f = 1.0;                                  // m!
        for (int m = 0; m <= 4; ++m) {
            const double sgn = (m % 2 == 0) ? 1.0 : -1.0;
            M[m + 1] = (sgn * f / std::pow(om, m + 1.0) - (m + 1.0) * M[m]) / th;
            f *= (m + 1.0);
        }
    }
}

static inline void d5nb1_row(double y, double mu, double th, double* out) {
    if (th < 1.0) {
        // the form in mu + i theta
        double S[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        for (double i = 0.0; i < y; i += 1.0) {
            const double t = mu + i * th;
            const double t2 = t * t;
            const double inv = 1.0 / (t2 * t2 * t);     // 1 / t^5
            double ib = 1.0;                             // i^b
            for (int b = 0; b <= 5; ++b) {
                S[b] += ib * inv;
                ib *= i;
            }
        }
        double M[6];
        d5nb1_M(th, M);
        const double om = 1.0 + th;
        const double om5 = om * om * om * om * om;
        for (int k = 0; k <= 5; ++k) out[k] = 24.0 * S[k];
        out[4] -= M[4];                                  // one mu: -M_4
        out[5] -= mu * M[5] + 24.0 * y / om5;            // no mu: -mu M_5 - 24 y/(1+theta)^5
        return;
    }

    // the size route
    const double r = mu / th;
    const double D0 = d5nb1_psi_shift(0, y, r);
    const double D1 = d5nb1_psi_shift(1, y, r);
    const double D2 = d5nb1_psi_shift(2, y, r);
    const double D3 = d5nb1_psi_shift(3, y, r);
    const double D4 = d5nb1_psi_shift(4, y, r);

    const double Lg = std::log1p(th);
    const double th2 = th * th, th3 = th2 * th, th4 = th2 * th2;
    const double th5 = th4 * th, th6 = th5 * th, th7 = th6 * th, th8 = th7 * th;
    const double th9 = th8 * th, th10 = th9 * th;
    const double om = 1.0 + th;
    const double om4 = om * om * om * om;
    const double om5 = om4 * om;
    const double mu2 = mu * mu, mu3 = mu2 * mu, mu4 = mu3 * mu, mu5 = mu4 * mu;

    // elementary parts, from derive.py, with log(1+theta) as Lg
    const double elem4 = -2.0 * (12.0 * th4 * Lg - 25.0 * th4 + 48.0 * th3 * Lg
                                 - 52.0 * th3 + 72.0 * th2 * Lg - 42.0 * th2
                                 + 48.0 * th * Lg - 12.0 * th + 12.0 * Lg)
                         / (th5 * om4);
    const double elem5 = 2.0 * (60.0 * mu * th5 * Lg - 137.0 * mu * th5
                                + 300.0 * mu * th4 * Lg - 385.0 * mu * th4
                                + 600.0 * mu * th3 * Lg - 470.0 * mu * th3
                                + 600.0 * mu * th2 * Lg - 270.0 * mu * th2
                                + 300.0 * mu * th * Lg - 60.0 * mu * th
                                + 60.0 * mu * Lg
                                + 60.0 * th5 * y + 120.0 * th4 * y + 120.0 * th3 * y
                                + 60.0 * th2 * y + 12.0 * th * y)
                         / (th6 * om5);

    out[0] = D4 / th5;
    out[1] = -4.0 * D3 / th5 - mu * D4 / th6;
    out[2] = 12.0 * D2 / th5 + 8.0 * mu * D3 / th6 + mu2 * D4 / th7;
    out[3] = -24.0 * D1 / th5 - 36.0 * mu * D2 / th6 - 12.0 * mu2 * D3 / th7
             - mu3 * D4 / th8;
    out[4] = elem4 + 24.0 * D0 / th5 + 96.0 * mu * D1 / th6 + 72.0 * mu2 * D2 / th7
             + 16.0 * mu3 * D3 / th8 + mu4 * D4 / th9;
    out[5] = elem5 - 120.0 * mu * D0 / th6 - 240.0 * mu2 * D1 / th7
             - 120.0 * mu3 * D2 / th8 - 20.0 * mu4 * D3 / th9 - mu5 * D4 / th10;
}

// [[Rcpp::export]]
List negbin1_deriv5_cpp(NumericVector y, NumericVector mu, NumericVector theta) {
    const int ny = y.size(), nm = mu.size(), nt = theta.size();
    const int n = std::max(ny, std::max(nm, nt));
    NumericVector m_mmmmm(n), m_mmmmt(n), m_mmmtt(n), m_mmttt(n), m_mtttt(n), m_tttttt(n);

    for (int i = 0; i < n; ++i) {
        double out[6];
        d5nb1_row(y[i % ny], mu[i % nm], theta[i % nt], out);
        m_mmmmm[i] = out[0];
        m_mmmmt[i] = out[1];
        m_mmmtt[i] = out[2];
        m_mmttt[i] = out[3];
        m_mtttt[i] = out[4];
        m_tttttt[i] = out[5];
    }

    return List::create(
        Named("mu_mu_mu_mu_mu") = m_mmmmm,
        Named("mu_mu_mu_mu_theta") = m_mmmmt,
        Named("mu_mu_mu_theta_theta") = m_mmmtt,
        Named("mu_mu_theta_theta_theta") = m_mmttt,
        Named("mu_theta_theta_theta_theta") = m_mtttt,
        Named("theta_theta_theta_theta_theta") = m_tttttt
    );
}
