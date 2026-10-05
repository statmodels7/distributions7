#include <Rcpp.h>
#include "d7_par.h"
#include "psi_diff.h"
#include "pt_negbin1.h"
#include "pt_sqrt.h"
using namespace Rcpp;

// Negative binomial with a variance LINEAR in the mean: Var(Y) = mu (1 + theta),
// against the quadratic mu + mu^2/theta of negbin_distrib(). The two are
// different families rather than two parametrizations of one, and the
// difference shows in where the mean sits: here the size is r = mu/theta, so
// mu appears INSIDE the gamma functions, while in NB2 it stays outside them.
//
// With r = mu/theta and p = 1/(1+theta),
//   l = lgamma(y+r) - lgamma(r) - lgamma(y+1) - r log(1+theta)
//       + y log(theta) - y log(1+theta).
// Writing P = dl/dr and Q = dl/dtheta at fixed r,
//   P    = psi(y+r) - psi(r) - log(1+theta),
//   Q    = -r/(1+theta) + y/theta - y/(1+theta),
//   P_r  = psi'(y+r) - psi'(r),
//   P_th = -1/(1+theta)  (which is also Q_r, the mixed second derivative),
//   Q_th = r/(1+theta)^2 - y/theta^2 + y/(1+theta)^2,
// and the chain rule through r = mu/theta gives the rest.

// NB1parts, nb1_parts() and nb1_E_Pr() are in pt_negbin1.h.

// The three components of the expected hessian at one (mu, theta), written
// apart because the caller reads them from two places.
static inline void nb1_E_parts(double m, double t,
                               double &emm, double &emt, double &ett) {
    double om = 1.0 + t;
    double EPr = nb1_E_Pr(m, t);
    double rm = 1.0 / t, rt = -m / (t * t);
    double Pth = -1.0 / om;
    emm = d7::negbin1_expected_mu_mu(m, t, EPr);
    emt = EPr * rm * rt + Pth * rm;
    ett = d7::negbin1_expected_theta_theta(m, t, EPr);
    // r_tt appears nowhere: the term it multiplies carries P, whose mean is
    // zero by the first Bartlett identity.
}

// [[Rcpp::export]]
NumericVector negbin1_logpmf_cpp(NumericVector y, NumericVector mu,
                                 NumericVector theta,
                        int threads = 1) {
    int n = y.size();
    NumericVector out(n);
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_s ? mu[0] : mu[i];
        double t = th_s ? theta[0] : theta[i];
        out[i] = R::dnbinom(y[i], m / t, 1.0 / (1.0 + t), 1);
    });
    return out;
}

// [[Rcpp::export]]
List negbin1_gradient_cpp(NumericVector y, NumericVector mu, NumericVector theta,
                        int threads = 1) {
    int n = y.size();
    NumericVector g_mu(n), g_th(n);
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_s ? mu[0] : mu[i];
        double t = th_s ? theta[0] : theta[i];
        NB1parts z = nb1_parts(y[i], m, t, 1);
        g_mu[i] = d7::negbin1_score_mu(z);
        g_th[i] = d7::negbin1_score_theta(z);
    });
    return List::create(Named("mu") = g_mu, Named("theta") = g_th);
}

// [[Rcpp::export]]
List negbin1_hessian_cpp(NumericVector y, NumericVector mu, NumericVector theta,
                        int threads = 1) {
    int n = y.size();
    NumericVector h_mm(n), h_mt(n), h_tt(n);
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_s ? mu[0] : mu[i];
        double t = th_s ? theta[0] : theta[i];
        NB1parts z = nb1_parts(y[i], m, t, 2);
        h_mm[i] = d7::negbin1_hess_mu_mu(z);
        h_mt[i] = (z.Pr * z.rt + z.Pth) * z.rm + z.P * z.rmt;
        h_tt[i] = d7::negbin1_hess_theta_theta(z);
    });
    return List::create(Named("mu_mu") = h_mm,
                        Named("mu_theta") = h_mt,
                        Named("theta_theta") = h_tt);
}

// The expected Hessian. E[P] = 0 by the first Bartlett identity -- the score
// in mu is P/theta, so its mean vanishing means P's does -- which removes
// every term carrying P and leaves E[psi'(Y+r)] as the only quantity without
// a closed form.
//
// One observation is computed and written in full by one thread, as every
// other kernel here does. What this replaced was a sequential loop that
// memoized the expectation across CONSECUTIVE equal parameters: that fires
// wherever a design repeats a mean in adjacent rows, and nowhere at all
// under an offset, which gives every row its own -- measured on the same
// 28800 cells, 28800 distinct means of 28800, so the memo never hit once
// and the loop could not be parallelized because of it.
// [[Rcpp::export]]
List negbin1_expected_hessian_cpp(NumericVector y, NumericVector mu,
                                  NumericVector theta,
                        int threads = 1) {
    int n = y.size();
    NumericVector h_mm(n), h_mt(n), h_tt(n);
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);
    bool both_scalar = mu_s && th_s;

    double emm0 = 0, emt0 = 0, ett0 = 0;
    if (both_scalar) nb1_E_parts(mu[0], theta[0], emm0, emt0, ett0);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double emm = emm0, emt = emt0, ett = ett0;
        if (!both_scalar) {
            double m = mu_s ? mu[0] : mu[i];
            double t = th_s ? theta[0] : theta[i];
            nb1_E_parts(m, t, emm, emt, ett);
        }
        h_mm[i] = emm; h_mt[i] = emt; h_tt[i] = ett;
    });
    return List::create(Named("mu_mu") = h_mm,
                        Named("mu_theta") = h_mt,
                        Named("theta_theta") = h_tt);
}

// ---- the derivatives of the expected information ------------------------
//
// Every component is c0(mu, theta) + c1(mu, theta) G with c0, c1 rational
// (nb1_E_parts above) and G = E[psi'(Y + r) - psi'(r)], r = mu/theta, the one
// quantity without a closed form. G is a sum over the support whose mass moves
// with the parameters, so in the coordinates (r, theta) its derivatives are
// sums over the same mass with the score of the mass beside the summand:
// with a = d log p/dr = S - log1p(theta), S = sum_{j<y} 1/(r+j),
// b = d log p/dtheta = (y - mu)/(theta (1 + theta)) and
// b' = d2 log p/dtheta2 = (r + y)/(1 + theta)^2 - y/theta^2,
//   G      = E[T]
//   G_r    = E[a T + U]
//   G_t    = E[b T]
//   G_rr   = E[a^2 T + T^2 + 2 a U + V]
//   G_rt   = E[b (a T + U) - T/(1 + theta)]
//   G_tt   = E[(b^2 + b') T]
// where T, U, V are the differences psi^(n)(r + y) - psi^(n)(r) for n = 1, 2,
// 3 by their exact recurrences, T = -sum 1/(r+j)^2, U = sum 2/(r+j)^3,
// V = -sum 6/(r+j)^4, accumulated beside the mass as nb1_E_Pr accumulates T.
// The mass is the one nb1_E_Pr sums, seeded and switched out of log scale by
// the same rule, and the loop stops by the same rule on its terms.
// The first order reads G, G_r and G_t, which need T and U beside the mass;
// the second all six, which need V as well.
// nb1_G_derivs1() is in pt_negbin1.h; the second order follows it.
static void nb1_G_derivs2(double mu, double th, double* G) {
    double r = mu / th;
    double ratio = th / (1.0 + th);
    double lratio = std::log(th) - std::log1p(th);
    double L = std::log1p(th), op = 1.0 + th, iop2 = 1.0 / (op * op);
    double cap = 100.0 + mu + 20.0 * d7::sqrt_cr(mu * (1.0 + th))
                 + 80.0 / (-std::log(ratio));
    int kmax = (int) std::min(cap, 2.0e9);
    const double mode = negbin1_series_mode(r, th);
    double lpk = -r * std::log1p(th);
    bool logscale = (lpk <= -640.0);
    double pk = std::exp(lpk);
    double S = 0.0, T = 0.0, U = 0.0, V = 0.0, cum = 0.0;
    double s[6] = {0, 0, 0, 0, 0, 0};
    for (int k = 0; k <= kmax; ++k) {
        double kd = (double) k;
        double a = S - L;
        double b = (kd - mu) / (th * op);
        double bt = (r + kd) * iop2 - kd / (th * th);
        s[0] += pk * T;
        s[1] += pk * (a * T + U);
        s[2] += pk * b * T;
        s[3] += pk * (a * a * T + T * T + 2.0 * a * U + V);
        s[4] += pk * (b * (a * T + U) - T / op);
        s[5] += pk * (b * b + bt) * T;
        cum += pk;
        if (negbin1_series_done(kd, pk, mode)) break;
        double v = 1.0 / (r + kd), v2 = v * v;
        S += v;
        T -= v2;
        U += 2.0 * v2 * v;
        V -= 6.0 * v2 * v2;
        if (logscale) {
            lpk += lratio + std::log((kd + r) / (kd + 1.0));
            pk = std::exp(lpk);
            if (lpk > -640.0) logscale = false;
        } else {
            pk *= (kd + r) / (kd + 1.0) * ratio;
        }
    }
    for (int j = 0; j < 6; ++j) G[j] = (cum > 0) ? s[j] / cum : 0.0;
}

// [[Rcpp::export]]
List negbin1_dexpected1_cpp(NumericVector y, NumericVector mu,
                            NumericVector theta, int threads = 1) {
    int n = y.size();
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);
    const int K = 6;
    std::vector<NumericVector> v(K);
    std::vector<double*> p(K);
    for (int k = 0; k < K; ++k) { v[k] = NumericVector(n); p[k] = v[k].begin(); }
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_s ? mu[0] : mu[i];
        double t = th_s ? theta[0] : theta[i];
        double r = m / t, g[3];
        nb1_G_derivs1(m, t, g);
        // G and its first derivatives in (mu, theta): d_mu = d_r / t,
        // d_theta|mu = d_theta|r - (r/t) d_r
        double G = g[0], it = 1.0 / t, rt = r * it;
        double Gm = g[1] * it, Gt = g[2] - rt * g[1];
        double it2 = it * it, it3 = it2 * it, it4 = it2 * it2, it5 = it4 * it;
        double op = 1.0 + t, iop = 1.0 / op, iop2 = iop * iop, iop3 = iop2 * iop;
        // c0, c1 and their first derivatives, per component (mm, tt, mt):
        // index 0 value, 1 d_mu, 2 d_theta
        // E_mm = G / t^2
        double a0[3] = {it2, 0.0, -2.0 * it3};
        double z0[3] = {0, 0, 0};
        // E_tt = m^2 G / t^4 + m h(t)
        double a1[3] = {m * m * it4, 2.0 * m * it4, -4.0 * m * m * it5};
        double h = 2.0 * it2 * iop + it * iop2 - it2 + iop2;
        double h1 = 2.0 * (-2.0 * it3 * iop - it2 * iop2)
                    - it2 * iop2 - 2.0 * it * iop3 + 2.0 * it3 - 2.0 * iop3;
        double z1[3] = {m * h, h, m * h1};
        // E_mt = -m G / t^3 - 1/(t (1 + t))
        double a2[3] = {-m * it3, -it3, 3.0 * m * it4};
        double q = t * op, q1 = 1.0 + 2.0 * t, iq = 1.0 / q, iq2 = iq * iq;
        double z2[3] = {-iq, 0.0, q1 * iq2};
        const double* A[3] = {a0, a1, a2};
        const double* Z[3] = {z0, z1, z2};
        // the two diagonal components, d_mu E_mm and d_theta E_tt, come
        // from pt_negbin1.h, which the registry reads
        for (int e = 0; e < 3; ++e) {
            if (e != 0) p[2 * e][i] = Z[e][1] + A[e][1] * G + A[e][0] * Gm;
            if (e != 1) p[2 * e + 1][i] = Z[e][2] + A[e][2] * G + A[e][0] * Gt;
        }
        p[0][i] = d7::negbin1_dexpected_mu_mu_mu(m, t, g);
        p[3][i] = d7::negbin1_dexpected_theta_theta_theta(m, t, g);
    });
    CharacterVector nm = CharacterVector::create(
        "mu_mu_mu", "mu_mu_theta", "theta_theta_mu",
        "theta_theta_theta", "mu_theta_mu", "mu_theta_theta");
    List out(K);
    for (int k = 0; k < K; ++k) out[k] = v[k];
    out.attr("names") = nm;
    return out;
}

// [[Rcpp::export]]
List negbin1_dexpected2_cpp(NumericVector y, NumericVector mu,
                            NumericVector theta, int threads = 1) {
    int n = y.size();
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);
    const int K = 9;
    std::vector<NumericVector> v(K);
    std::vector<double*> p(K);
    for (int k = 0; k < K; ++k) { v[k] = NumericVector(n); p[k] = v[k].begin(); }
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_s ? mu[0] : mu[i];
        double t = th_s ? theta[0] : theta[i];
        double r = m / t, g[6];
        nb1_G_derivs2(m, t, g);
        // G and its derivatives in (mu, theta): d_mu = d_r / t,
        // d_theta|mu = d_theta|r - (r/t) d_r
        double G = g[0], it = 1.0 / t, rt = r * it;
        double Gm = g[1] * it, Gt = g[2] - rt * g[1];
        double Gmm = g[3] * it * it;
        double Gmt = (g[4] - rt * g[3]) * it - g[1] * it * it;
        double Gtt = g[5] - 2.0 * rt * g[4] + rt * rt * g[3] + 2.0 * r * it * it * g[1];
        double it2 = it * it, it3 = it2 * it, it4 = it2 * it2, it5 = it4 * it,
               it6 = it4 * it2;
        double op = 1.0 + t, iop = 1.0 / op, iop2 = iop * iop, iop3 = iop2 * iop,
               iop4 = iop2 * iop2;
        // c0, c1 and their derivatives, per component (mm, tt, mt):
        // index 0 value, 1 d_mu, 2 d_theta, 3 mu mu, 4 theta theta, 5 mu theta
        double c1[3][6], c0[3][6];
        // E_mm = G / t^2
        double a0[6] = {it2, 0.0, -2.0 * it3, 0.0, 6.0 * it4, 0.0};
        double z0[6] = {0, 0, 0, 0, 0, 0};
        // E_tt = m^2 G / t^4 + m h(t)
        double a1[6] = {m * m * it4, 2.0 * m * it4, -4.0 * m * m * it5,
                        2.0 * it4, 20.0 * m * m * it6, -8.0 * m * it5};
        double h = 2.0 * it2 * iop + it * iop2 - it2 + iop2;
        double h1 = 2.0 * (-2.0 * it3 * iop - it2 * iop2)
                    - it2 * iop2 - 2.0 * it * iop3 + 2.0 * it3 - 2.0 * iop3;
        double h2 = 2.0 * (6.0 * it4 * iop + 4.0 * it3 * iop2 + 2.0 * it2 * iop3)
                    + 2.0 * it3 * iop2 + 4.0 * it2 * iop3 + 6.0 * it * iop4
                    - 6.0 * it4 + 6.0 * iop4;
        double z1[6] = {m * h, h, m * h1, 0.0, m * h2, h1};
        // E_mt = -m G / t^3 - 1/(t (1 + t))
        double a2[6] = {-m * it3, -it3, 3.0 * m * it4, 0.0, -12.0 * m * it5, 3.0 * it4};
        double q = t * op, q1 = 1.0 + 2.0 * t, iq = 1.0 / q, iq2 = iq * iq;
        double z2[6] = {-iq, 0.0, q1 * iq2, 0.0, 2.0 * iq2 - 2.0 * q1 * q1 * iq2 * iq, 0.0};
        for (int j = 0; j < 6; ++j) {
            c1[0][j] = a0[j]; c0[0][j] = z0[j];
            c1[1][j] = a1[j]; c0[1][j] = z1[j];
            c1[2][j] = a2[j]; c0[2][j] = z2[j];
        }
        for (int e = 0; e < 3; ++e) {
            double *A = c1[e], *Z = c0[e];
            p[3 * e][i] = Z[3] + A[3] * G + 2.0 * A[1] * Gm + A[0] * Gmm;
            p[3 * e + 1][i] = Z[4] + A[4] * G + 2.0 * A[2] * Gt + A[0] * Gtt;
            p[3 * e + 2][i] = Z[5] + A[5] * G + A[1] * Gt + A[2] * Gm + A[0] * Gmt;
        }
    });
    CharacterVector nm = CharacterVector::create(
        "mu_mu_mu_mu", "mu_mu_theta_theta", "mu_mu_mu_theta",
        "theta_theta_mu_mu", "theta_theta_theta_theta",
        "theta_theta_mu_theta", "mu_theta_mu_mu",
        "mu_theta_theta_theta", "mu_theta_mu_theta");
    List out(K);
    for (int k = 0; k < K; ++k) out[k] = v[k];
    out.attr("names") = nm;
    return out;
}
