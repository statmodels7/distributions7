// The third and fourth derivatives of the expected information of the
// Poisson-inverse Gaussian, family pig2 (parameters mu, alpha), on the parameter
// scale. Three pieces, one per pair of the components.
//
// The pairs mu_mu and mu_alpha are closed forms. In this parametrization
// E_mu_mu = -1/mu + 1/sqrt(mu^2 + alpha^2) and E_mu_alpha = 0 (orthogonality) for
// every (mu, alpha); the identity was found in the series of pig_asym.h (the
// coefficients of pig2_asym0_0 are those of 1/sqrt(1 + t^2) in t = mu/alpha) and
// checked against the 50-digit sums of the reference at 18 points of mu in
// [0.05, 1e4], alpha in [0.05, 1e5], to 1e-24 of the order's largest component,
// for all components of both orders. The derivatives of the sum over the support
// of these two pairs cannot be evaluated accurately near the Poisson limit: the
// terms are of size (sd/mu)^k and the results are mu^-(k-1) times a function of
// mu/alpha, so that round-off of the terms alone is 1e-7 (order 3) and 1.5e-5
// (order 4) of the largest component at (mu, alpha) = (1e4, 1e5), and the mixed
// pair, which is zero, carries that noise as its value.
//
// The pair alpha_alpha has three routes.
//
// (1) The series. E_aa(mu, alpha) = sum_k sum_m c[k][m] mu^m alpha^-k, 4 <= k,
// 2 <= m <= k - 2 (the series of pig_asym.h, pig2_asym0_1, there to k = 25),
// extended to k = 70 in pig2_aa_series.h. Its derivatives of orders 3 and 4 are
// taken term by term, aa_series<K>() below. The series converges for t =
// (1 + mu)/alpha up to about 0.5 once mu is of order 100 or more (the ratio of
// successive terms is about t), and is asymptotic for small mu, where the minimal
// term grows with t. The sum is taken by diagonals k, which gives the truncation
// error: the sum is accepted when the last two diagonals are below 2e-11 of the
// block's largest component and decreasing, or when three consecutive diagonals
// are below 1e-17 of it. Otherwise the walk is used. The series needs no row of
// the log-density, so it has no cancellation: its error is 1e-15 of the block's
// largest component where it is accepted, against up to 1e-3 for the walk at
// mu = 1e3 near the Poisson limit.
//
// (1b) The same series summed by powers of 1/alpha at fixed u = mu/alpha (pig2_aa_rows.h).
// Each coefficient function S_n(u) of alpha^-n is algebraic, a polynomial over a power of
// 1 + u^2 for the even and for the odd powers of u (identified from the series for
// n = 2..23), so that the rows hold for any u, where the diagonal series needs t below
// 0.5. The rows are the terms of an expansion in 1/mu and 1/alpha: they are accepted once
// alpha is above about 150 (mu >= 300), 100 (mu >= 100) or 70 (mu >= 30), to 1e-15 of
// the block for mu >= 1e3 and alpha >= 500, whatever t. Tried when (1) is not accepted.
//
// (2) The sum over the support, used where neither is accepted (small mu with t above
// about 0.1, alpha below about 100 to 150). There its error is below 1e-14 of the block
// wherever it was measured (mu = 100, alpha = 10 to 100; mu = 1e3, alpha = 100).
// E_aa = E[l_aa] = -E[l_a l_a], and the derivatives are taken under the sum:
// d_c E[h] = E[D_c h] with D_c h = d_c h + h l_c, so that d_S E_aa =
// E[ D_S(-l_a l_a) ] for a multiset S of coordinates. Each component is a
// polynomial in the partials L{i}{j} = d^(i+j) l / dmu^i dalpha^j, derived
// symbolically (misure_d34e/agents/pig2/derive.py, accumulation in gen_aa.py);
// L{i}{j} with i + j <= 4 is read from the row of pt_pig_d34e.h and i + j = 5 from
// its fifth-order row. The row and the walk over the support are those of
// pig_expected_obs() in pig_hd.cpp, one walk per observation for all components
// of the order.
//
// The components are in the order of d3expected_names() and d4expected_names():
// the pairs (mu_mu, alpha_alpha, mu_alpha) in turn, and within a pair the index
// multisets with j indices on alpha, j = 0, 1, ...

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "d7_par.h"
#include "d34e_run.h"
#include "pt_pig_d34e.h"
#include "pt_pig2_aa.h"
using namespace Rcpp;

namespace {

using namespace pig2_aa;   // aa_series, aa_rows and pig2_mm_orders (pt_pig2_aa.h)

struct Acc3 {
    // the components of the alpha_alpha pair, order 3, for one observation,
    // accumulated into o[4 + j] (j indices on alpha); f is the mass
    void operator()(double f, const double* v, double* o) const {
    const double L01 = v[2];
    const double L02 = v[5];
    const double L03 = v[9];
    const double L04 = v[14];
    const double L10 = v[1];
    const double L11 = v[4];
    const double L12 = v[8];
    const double L13 = v[13];
    const double L20 = v[3];
    const double L21 = v[7];
    const double L22 = v[12];
    const double L30 = v[6];
    const double L31 = v[11];
    const double t0 = (2.0*L01);
    const double t1 = (6.0*L21);
    const double t2 = (L01*L01);
    const double t3 = (L01*L10);
    const double t4 = (6.0*L11);
    const double t5 = (L01*t4);
    const double t6 = (6.0*(L11*L11));
    const double t7 = (L10*L10);
    const double t8 = (3.0*t2);
    const double t9 = (4.0*L12);
    const double t10 = (L01*L01*L01);
    const double t11 = (L02*t0);
    const double t12 = (4.0*L11);
    const double t13 = (L02*L10);
    const double t14 = (5.0*t2);
    const double t15 = (L02*L02);
    o[4] += f * ((-1.0*L10*t6) + (-1.0*L11*t1) + (-1.0*L20*t5) + (-1.0*L30*t2) + (-1.0*L31*t0) + (-1.0*t1*t3) + (-1.0*t2*(L10*L10*L10)) + (-1.0*t5*t7) + (-1.0*L10*L20*t8));
    o[5] += f * ((-1.0*L01*t6) + (-1.0*L11*t9) + (-1.0*L20*t10) + (-1.0*L20*t11) + (-1.0*L21*t8) + (-1.0*L22*t0) + (-1.0*t10*t7) + (-1.0*t11*t7) + (-1.0*t12*t13) + (-1.0*t3*t9) + (-2.0*L02*L21) + (-1.0*L10*t2*t4));
    o[6] += f * ((-1.0*L02*t9) + (-1.0*L10*(L01*L01*L01*L01)) + (-1.0*L12*t14) + (-1.0*L13*t0) + (-1.0*t10*t12) + (-1.0*t13*t14) + (-2.0*L03*L11) + (-2.0*L10*t15) + (-1.0*L03*L10*t0) + (-10.0*L01*L02*L11));
    o[7] += f * ((-1.0*(L01*L01*L01*L01*L01)) + (-1.0*L04*t0) + (-12.0*L01*t15) + (-9.0*L02*t10) + (-7.0*L03*t2) + (-6.0*L02*L03));
    }
};


struct Acc4 {
    // the components of the alpha_alpha pair, order 4, for one observation,
    // accumulated into o[5 + j] (j indices on alpha); f is the mass
    void operator()(double f, const double* v, double* o) const {
    const double L01 = v[2];
    const double L02 = v[5];
    const double L03 = v[9];
    const double L04 = v[14];
    const double L05 = v[20];
    const double L10 = v[1];
    const double L11 = v[4];
    const double L12 = v[8];
    const double L13 = v[13];
    const double L14 = v[19];
    const double L20 = v[3];
    const double L21 = v[7];
    const double L22 = v[12];
    const double L23 = v[18];
    const double L30 = v[6];
    const double L31 = v[11];
    const double L32 = v[17];
    const double L40 = v[10];
    const double L41 = v[16];
    const double t0 = (2.0*L01);
    const double t1 = (8.0*L31);
    const double t2 = (L01*L01);
    const double t3 = (L01*L10);
    const double t4 = (L01*L11);
    const double t5 = (8.0*t4);
    const double t6 = (L01*L20);
    const double t7 = (12.0*L21);
    const double t8 = (L11*L21);
    const double t9 = (24.0*L10);
    const double t10 = (L11*L11);
    const double t11 = (12.0*t10);
    const double t12 = (L10*L10);
    const double t13 = (L01*t12);
    const double t14 = (L10*L10*L10);
    const double t15 = (L10*t2);
    const double t16 = (3.0*t2);
    const double t17 = (6.0*L22);
    const double t18 = (6.0*L12);
    const double t19 = (L01*L01*L01);
    const double t20 = (L02*t0);
    const double t21 = (6.0*L02);
    const double t22 = (L10*L21);
    const double t23 = (L11*t21);
    const double t24 = (12.0*L11);
    const double t25 = (L10*L12);
    const double t26 = (9.0*t2);
    const double t27 = (L11*t26);
    const double t28 = (L10*t19);
    const double t29 = (4.0*L13);
    const double t30 = (L01*L01*L01*L01);
    const double t31 = (10.0*L02);
    const double t32 = (L03*t0);
    const double t33 = (8.0*L02);
    const double t34 = (L03*L10);
    const double t35 = (5.0*t2);
    const double t36 = (L02*L02);
    const double t37 = (2.0*t36);
    const double t38 = (L02*L11);
    const double t39 = (L02*t35);
    const double t40 = (L01*L02);
    const double t41 = (7.0*t2);
    const double t42 = (9.0*t19);
    o[5] += f * ((-6.0*(L21*L21)) + (-1.0*L11*t1) + (-1.0*L20*t11) + (-1.0*L30*t5) + (-1.0*L40*t2) + (-1.0*L41*t0) + (-1.0*t1*t3) + (-1.0*t11*t12) + (-1.0*t13*t7) + (-1.0*t14*t5) + (-1.0*t16*(L20*L20)) + (-1.0*t2*(L10*L10*L10*L10)) + (-1.0*t6*t7) + (-1.0*t8*t9) + (-4.0*L30*t15) + (-1.0*L11*t6*t9) + (-6.0*L20*t12*t2));
    o[6] += f * ((-6.0*(L11*L11*L11)) + (-1.0*L11*t17) + (-1.0*L20*t23) + (-1.0*L20*t27) + (-1.0*L21*t18) + (-1.0*L30*t19) + (-1.0*L30*t20) + (-1.0*L31*t16) + (-1.0*L32*t0) + (-1.0*t12*t23) + (-1.0*t12*t27) + (-1.0*t13*t18) + (-1.0*t14*t19) + (-1.0*t14*t20) + (-1.0*t17*t3) + (-1.0*t18*t6) + (-1.0*t21*t22) + (-1.0*t22*t26) + (-1.0*t24*t25) + (-18.0*L01*t8) + (-18.0*t10*t3) + (-3.0*L20*t28) + (-2.0*L02*L31) + (-1.0*L10*t21*t6));
    o[7] += f * ((-4.0*(L12*L12)) + (-1.0*L11*t29) + (-1.0*L20*t30) + (-1.0*L20*t32) + (-1.0*L20*t37) + (-1.0*L20*t39) + (-1.0*L22*t35) + (-1.0*L23*t0) + (-1.0*t10*t31) + (-1.0*t11*t2) + (-1.0*t12*t30) + (-1.0*t12*t32) + (-1.0*t12*t37) + (-1.0*t12*t39) + (-1.0*t25*t33) + (-1.0*t29*t3) + (-20.0*L12*t4) + (-20.0*t3*t38) + (-10.0*L12*t15) + (-8.0*L11*t28) + (-4.0*L02*L22) + (-4.0*L11*t34) + (-4.0*L21*t19) + (-2.0*L03*L21) + (-1.0*L01*L21*t31));
    o[8] += f * ((-1.0*L03*t18) + (-1.0*L10*(L01*L01*L01*L01*L01)) + (-1.0*L12*t42) + (-1.0*L13*t21) + (-1.0*L13*t41) + (-1.0*L14*t0) + (-1.0*t21*t34) + (-1.0*t24*t36) + (-1.0*t34*t41) + (-27.0*t2*t38) + (-24.0*L12*t40) + (-14.0*L03*t4) + (-12.0*t3*t36) + (-5.0*L11*t30) + (-2.0*L04*L11) + (-1.0*L02*L10*t42) + (-1.0*L04*L10*t0));
    o[9] += f * ((-1.0*(L01*L01*L01*L01*L01*L01)) + (-12.0*(L02*L02*L02)) + (-6.0*(L03*L03)) + (-1.0*L04*t26) + (-1.0*L04*t33) + (-1.0*L05*t0) + (-44.0*L03*t40) + (-39.0*t2*t36) + (-16.0*L03*t19) + (-14.0*L02*t30));
    }
};

}  // namespace

// [[Rcpp::export]]
List pig2_dexpected3_cpp(NumericVector y, NumericVector mu, NumericVector alpha, int threads = 1) {
    const int n = (int) std::max(y.size(), std::max(mu.size(), alpha.size()));
    const bool one = mu.size() == 1 && alpha.size() == 1;
    const int nrun = one ? std::min(n, 1) : n;
    d34e::Par M(mu), X(alpha);
    List out = d34e::run(nrun, threads, d7::kMinCostly, 12,
                         [&](std::size_t i, d34e::Row o) {
        const double m = M[i], x = X[i];
        double c[12];
        bool ok = R_finite(m) && R_finite(x) && m > 0.0 && x > 0.0;
        if (ok) {
            const double t = m / x;
            const double b = x / (t + std::hypot(1.0, t));   // 1/sigma
            const double arg = x;
            const double tail = 2.0 * m / b;
            for (int k = 0; k < 12; ++k) c[k] = 0.0;
            const bool ser = aa_series<3>(m, x, c + 4) || aa_rows<3>(m, x, c + 4);
            ok = ser || d34e_pig::d34_walk<3, d7::PigSeq>(
                m, arg, tail,
                [&](double yy, const double* pre, double* v) {
                    d34e_pig::pig2_row(yy, m, x, 4, false, v, pre);
                },
                [&](double yy, double* v) {
                    d34e_pig::pig2_d5_row(yy, m, x, v + 15);
                },
                Acc3(), c);
            if (ok) {
                pig2_mm_orders(3, m, x, c);
                for (int k = 2 * 4; k < 12; ++k) c[k] = 0.0;   // E_mu_alpha = 0
            }
        }
        for (int k = 0; k < 12; ++k) o[k] = ok ? c[k] : NA_REAL;
    });
    if (nrun < n) {
        for (int k = 0; k < 12; ++k) {
            NumericVector col = out[k];
            out[k] = NumericVector(n, col[0]);
        }
    }
    CharacterVector nms(12);
    const char* nm[12] = {
        "mu_mu_mu_mu_mu", "mu_mu_mu_mu_alpha", "mu_mu_mu_alpha_alpha", "mu_mu_alpha_alpha_alpha",
        "alpha_alpha_mu_mu_mu", "alpha_alpha_mu_mu_alpha", "alpha_alpha_mu_alpha_alpha", "alpha_alpha_alpha_alpha_alpha",
        "mu_alpha_mu_mu_mu", "mu_alpha_mu_mu_alpha", "mu_alpha_mu_alpha_alpha", "mu_alpha_alpha_alpha_alpha",
    };
    for (int k = 0; k < 12; ++k) nms[k] = nm[k];
    out.names() = nms;
    return out;
}

// [[Rcpp::export]]
List pig2_dexpected4_cpp(NumericVector y, NumericVector mu, NumericVector alpha, int threads = 1) {
    const int n = (int) std::max(y.size(), std::max(mu.size(), alpha.size()));
    const bool one = mu.size() == 1 && alpha.size() == 1;
    const int nrun = one ? std::min(n, 1) : n;
    d34e::Par M(mu), X(alpha);
    List out = d34e::run(nrun, threads, d7::kMinCostly, 15,
                         [&](std::size_t i, d34e::Row o) {
        const double m = M[i], x = X[i];
        double c[15];
        bool ok = R_finite(m) && R_finite(x) && m > 0.0 && x > 0.0;
        if (ok) {
            const double t = m / x;
            const double b = x / (t + std::hypot(1.0, t));   // 1/sigma
            const double arg = x;
            const double tail = 2.0 * m / b;
            for (int k = 0; k < 15; ++k) c[k] = 0.0;
            const bool ser = aa_series<4>(m, x, c + 5) || aa_rows<4>(m, x, c + 5);
            ok = ser || d34e_pig::d34_walk<4, d7::PigSeq>(
                m, arg, tail,
                [&](double yy, const double* pre, double* v) {
                    d34e_pig::pig2_row(yy, m, x, 4, false, v, pre);
                },
                [&](double yy, double* v) {
                    d34e_pig::pig2_d5_row(yy, m, x, v + 15);
                },
                Acc4(), c);
            if (ok) {
                pig2_mm_orders(4, m, x, c);
                for (int k = 2 * 5; k < 15; ++k) c[k] = 0.0;   // E_mu_alpha = 0
            }
        }
        for (int k = 0; k < 15; ++k) o[k] = ok ? c[k] : NA_REAL;
    });
    if (nrun < n) {
        for (int k = 0; k < 15; ++k) {
            NumericVector col = out[k];
            out[k] = NumericVector(n, col[0]);
        }
    }
    CharacterVector nms(15);
    const char* nm[15] = {
        "mu_mu_mu_mu_mu_mu", "mu_mu_mu_mu_mu_alpha", "mu_mu_mu_mu_alpha_alpha", "mu_mu_mu_alpha_alpha_alpha",
        "mu_mu_alpha_alpha_alpha_alpha", "alpha_alpha_mu_mu_mu_mu", "alpha_alpha_mu_mu_mu_alpha", "alpha_alpha_mu_mu_alpha_alpha",
        "alpha_alpha_mu_alpha_alpha_alpha", "alpha_alpha_alpha_alpha_alpha_alpha", "mu_alpha_mu_mu_mu_mu", "mu_alpha_mu_mu_mu_alpha",
        "mu_alpha_mu_mu_alpha_alpha", "mu_alpha_mu_alpha_alpha_alpha", "mu_alpha_alpha_alpha_alpha_alpha",
    };
    for (int k = 0; k < 15; ++k) nms[k] = nm[k];
    out.names() = nms;
    return out;
}
