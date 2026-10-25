// The third and fourth derivatives of the expected information of the elastic
// net in its location mu, rate lambda and mixing weight alpha, on the
// parameter scale, in closed form. Each order is its own kernel
// (enet_dexpected3_cpp, enet_dexpected4_cpp), and both read the same derivatives
// of log Z in (a, c) through d34e_enet_table(), which is the recurrence of
// pt_enet_logz.h carried to order six (that header stops at order four).
//
// With a = lambda alpha, c = lambda (1 - alpha), Z(a, c) the normalizing
// constant and U = 1/Z the density at zero, the expected information does not
// depend on mu, so every derivative in mu is zero, and E_ml = E_ma = 0. The
// remaining components are
//   E_mm = -c - 2 a U                               (the form of pt_enet.h)
//   E_ll = -(alpha^2 D20 + 2 alpha bl D11 + bl^2 D02),  bl = 1 - alpha
//   E_aa = -lambda^2 (D20 - 2 D11 + D02)
//   E_la = -lambda (alpha D20 + (1 - 2 alpha) D11 - bl D02)
// with D_ij = d^(i+j) log Z / da^i dc^j. Each derivative in (lambda, alpha)
// moves a D_ij to D_(i+1)j or D_i(j+1), and U to -U (log Z)'. The derivatives
// are written out symbolically (misure_d34e/agents/enet/derive.py, which
// collects every coefficient in (lambda, alpha) in factored form) and printed
// by gen_cpp.py. Derivatives of log Z of order three and four therefore need
// D up to order five and six (for orders three and four of E_ll), which is why
// the two kernels build tables to those orders.
//
// The components are in the order of d3expected_names() and d4expected_names()
// for enet: for each pair (mu_mu, lambda_lambda, alpha_alpha, mu_lambda,
// mu_alpha, lambda_alpha), the multisets of derivative indices in the order
// they are listed there.
//
// Cancellation. The derivatives of log Z in (a, c) are formed from the moments
// of the density proportional to exp(-x s - s^2/2) on s > 0 (x = a/sqrt(c)), as
// in pt_enet_logz.h: the cumulants of the log Mills ratio for x < 2, and the
// cumulants of W = -U^2/2 with t = c/a^2 for x >= 2. No difference of large
// numbers is formed in the tables. The expressions of the components are
// polynomials in (lambda, alpha) times the D_ij and U, with the coefficients
// factored, so that no expanded polynomial in alpha is evaluated.
//
// The pair mu_mu does not read the table. E_mm = -c - 2 a U is a function of
// (lambda, alpha) alone, and its derivatives through D_ij and U cancel
// catastrophically when c is small (errors of 4e-8 at order 4 for lambda =
// 0.001, alpha = 1 - 1e-6). d34e_enet_emm() writes E_mm = -a^2 - 2 c + c m2(x),
// with m2 the second moment of the density proportional to exp(-x s - s^2/2),
// and differentiates the remainder c m2 by Taylor composition for x < 15 and
// the whole of E_mm by its asymptotic series in c/a^2 for x >= 15.
//
// The (lambda, alpha) pairs are found once per call (d34e_enet_pairs()), so the
// tables are formed per distinct pair, as enet.cpp does for the first
// derivatives. A pair of scalars costs one table.

#include <Rcpp.h>
#include <Rmath.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <utility>
#include <vector>
#include "d7_par.h"
#include "d34e_run.h"
using namespace Rcpp;

namespace {

// m_0 .. m_n of the density proportional to exp(-x s - s^2/2) on s > 0, x > 0,
// as enet_s_moments() of pt_enet_logz.h, with room for n up to twelve
inline void d34e_enet_moments(double x, int n, double* m) {
    m[0] = 1.0;
    if (x < 1.0) {
        // M_0 = Phi(-x)/phi(x), and M_(k+1) = k M_(k-1) - x M_k forward
        const double M0 = std::exp(R::pnorm(-x, 0.0, 1.0, 1, 1) -
                                   R::dnorm(x, 0.0, 1.0, 1));
        double Mm = M0, Mk = 1.0 - x * M0;
        if (n >= 1) m[1] = Mk / M0;
        for (int k = 1; k < n; ++k) {
            const double Mn = k * Mm - x * Mk;
            Mm = Mk;
            Mk = Mn;
            m[k + 1] = Mk / M0;
        }
        return;
    }
    // the continued fraction r_k = k/(x + r_(k+1)), run down from K
    const int K = x < 2.0 ? 2000 : (x < 10.0 ? 300 : 80);
    double r = 0.0, rr[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int k = K; k >= 1; --k) {
        r = k / (x + r);
        if (k <= n) rr[k] = r;
    }
    for (int k = 1; k <= n; ++k) m[k] = m[k - 1] * rr[k];
}

// kappa_1 .. kappa_n from the raw moments mu_0 = 1 .. mu_n (enet_cumulants())
inline void d34e_enet_cumul(const double* mu, int n, double* kap) {
    for (int m = 1; m <= n; ++m) {
        double s = mu[m];
        double binom = 1.0;  // choose(m - 1, i - 1)
        for (int i = 1; i < m; ++i) {
            s -= binom * kap[i] * mu[m - i];
            binom = binom * (m - i) / i;
        }
        kap[m] = s;
    }
}

// one Euler step on b (length len): b'_m = -h (j + m) b_m - q b_m,
// b'_(m + 1) += -h b_m; returns the new length (enet_euler_step())
inline int d34e_enet_euler(double* b, int len, int j, double h, int q) {
    double nb[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int m = 0; m < len; ++m) {
        nb[m] += (-h * (j + m) - q) * b[m];
        nb[m + 1] += -h * b[m];
    }
    for (int m = 0; m <= len; ++m) b[m] = nb[m];
    return len + 1;
}

// D[i][j] = d^(i + j) log Z / da^i dc^j for 1 <= i + j <= N (N <= 6), and
// *U = 1/Z, the density at zero. The branches and the Euler contractions are
// those of enet_logz_ac() in pt_enet_logz.h; only the order is carried on.
inline void d34e_enet_table(double a, double c, int N, double D[7][7],
                            double* U) {
    for (int i = 0; i < 7; ++i)
        for (int j = 0; j < 7; ++j) D[i][j] = 0.0;
    const double x = a / std::sqrt(c);
    double m[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    if (x >= 2.0) {
        // log Z = log 2 - log a + F(t), t = c/a^2, F^(j) the cumulants of W
        const double t = 1.0 / (x * x);
        double mu[9] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double F[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        d34e_enet_moments(x, 2 * N, m);
        double h = 1.0, x2 = 1.0;
        for (int j = 1; j <= N; ++j) {
            h *= -0.5;
            x2 *= x * x;
            mu[j] = h * x2 * m[2 * j];
        }
        d34e_enet_cumul(mu, N, F);
        for (int ord = 1; ord <= N; ++ord) {
            for (int j = 0; j <= ord; ++j) {
                const int i = ord - j;
                double b[9] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                int len = 1;
                for (int q = 0; q < i; ++q)
                    len = d34e_enet_euler(b, len, j, 2.0, q);
                double v = 0.0, tm = 1.0;
                for (int mm = 0; mm < len; ++mm) {
                    if (b[mm] != 0.0) v += b[mm] * tm * F[j + mm];
                    tm *= t;
                }
                v *= std::pow(a, -2.0 * j - i);
                if (j == 0) {
                    double f = 1.0;
                    for (int k = 2; k < i; ++k) f *= k;
                    v += ((i % 2) ? -f : f) * std::pow(a, -i);
                }
                D[i][j] = v;
            }
        }
    } else {
        // log Z = log 2 + L(x) - log(c)/2, L = log M_0
        double L[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double kap[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        d34e_enet_moments(x, N, m);
        d34e_enet_cumul(m, N, kap);
        for (int k = 1; k <= N; ++k) L[k] = (k % 2) ? -kap[k] : kap[k];
        for (int ord = 1; ord <= N; ++ord) {
            for (int j = 0; j <= ord; ++j) {
                const int i = ord - j;
                double b[9] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                int len = 1;
                for (int q = 0; q < j; ++q)
                    len = d34e_enet_euler(b, len, i, 0.5, q);
                double v = 0.0, xm = 1.0;
                for (int mm = 0; mm < len; ++mm) {
                    if (i + mm >= 1 && b[mm] != 0.0) v += b[mm] * xm * L[i + mm];
                    xm *= x;
                }
                v *= std::pow(c, -0.5 * i - j);
                if (i == 0) {
                    double f = 1.0;
                    for (int k = 2; k < j; ++k) f *= k;
                    v += -0.5 * ((j % 2) ? f : -f) * std::pow(c, -j);
                }
                D[i][j] = v;
            }
        }
    }
    // 1/Z = (a + sqrt(c) m_1)/2 with m_1 = M_1/M_0 (enet_density0())
    double m1[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    d34e_enet_moments(x, 1, m1);
    *U = 0.5 * (a + std::sqrt(c) * m1[1]);
}

// The distinct (lambda, alpha) pairs of the first n rows, with the tables of
// each to order N, and the index of each row's pair. tabD holds 49 entries
// per pair, tabD[u * 49 + 7 i + j] = D_ij of pair u.
inline void d34e_enet_pairs(const NumericVector& lambda,
                            const NumericVector& alpha, std::size_t n, int N,
                            std::vector<double>& ul, std::vector<double>& ua,
                            std::vector<double>& tabD,
                            std::vector<double>& tabU, std::vector<int>& id) {
    const d34e::Par L(lambda), A(alpha);
    std::map<std::pair<double, double>, int> seen;
    ul.clear();
    ua.clear();
    id.assign(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        const double l = L[i], al = A[i];
        const std::pair<double, double> key(l, al);
        auto it = seen.find(key);
        if (it == seen.end()) {
            it = seen.emplace(key, static_cast<int>(ul.size())).first;
            ul.push_back(l);
            ua.push_back(al);
        }
        id[i] = it->second;
    }
    tabD.assign(ul.size() * 49, 0.0);
    tabU.assign(ul.size(), 0.0);
    for (std::size_t u = 0; u < ul.size(); ++u) {
        double D[7][7];
        double U = 0.0;
        d34e_enet_table(ul[u] * ua[u], ul[u] * (1.0 - ua[u]), N, D, &U);
        for (int i = 0; i < 7; ++i)
            for (int j = 0; j < 7; ++j) tabD[u * 49 + 7 * i + j] = D[i][j];
        tabU[u] = U;
    }
}

// The derivatives of order D (3 or 4) of E_mm = -c - 2 a U, in (lambda, alpha),
// for the pair mu_mu. With x = a/sqrt(c) and U = sqrt(c) h(x)/2, h(x) =
// phi(x)/Phi(-x), E_mm = -c (1 + x h(x)) = -a^2 - 2 c + c m2(x), where m2(x) is
// the second moment of the density proportional to exp(-x s - s^2/2) on s > 0
// (x m1 + m2 = 1). The part -a^2 - 2 c is differentiated exactly. The
// derivatives of the table D_ij are not used for this pair: they cancel
// catastrophically once c is small (the terms are of the size of 1/c^k and
// the entry is finite), which was 4e-8 at (lambda, alpha) = (0.001, 1 - 1e-6).
// Two routes, by x:
//   x < 15:  the remainder c m2(x(lambda, alpha)) is composed as a bivariate
//            Taylor polynomial of total degree D about (lambda, alpha), with
//            m2^(k)(x) from the moments m_1 .. m_6 (d/dx m_j = m_j m_1 -
//            m_(j+1), so m2^(k) is a polynomial in the moments);
//   x >= 15: the asymptotic series E_mm = -sum_n phi_n c^n a^(2 - 2n), t =
//            c/a^2, phi(t) = t + 1/S(t), S(t) = sum_n (-1)^n (2n - 1)!! t^n,
//            each monomial in (lambda, alpha) differentiated exactly.

// m2^(k)(x), k = 0..4, from the moments m[0..6]
inline void d34e_enet_m2der(const double* m, double* d) {
    const double m1 = m[1], m2 = m[2], m3 = m[3], m4 = m[4], m5 = m[5], m6 = m[6];
    d[0] = m2;
    d[1] = m1*m2 - m3;
    d[2] = 2*((m1)*(m1))*m2 - 2*m1*m3 - ((m2)*(m2)) + m4;
    d[3] = 6*((m1)*(m1)*(m1))*m2 - 6*((m1)*(m1))*m3 - 6*m1*((m2)*(m2)) + 3*m1*m4 + 4*m2*m3 - m5;
    d[4] = 24*((m1)*(m1)*(m1)*(m1))*m2 - 24*((m1)*(m1)*(m1))*m3 - 36*((m1)*(m1))*((m2)*(m2)) + 12*((m1)*(m1))*m4 + 32*m1*m2*m3 - 4*m1*m5 + 6*((m2)*(m2)*(m2)) - 7*m2*m4 - 4*((m3)*(m3)) + m6;
}

// falling factorial q (q - 1) ... (q - s + 1)
inline double d34e_enet_ff(double q, int s) {
    double r = 1.0;
    for (int k = 0; k < s; ++k) r *= (q - k);
    return r;
}

// C = A B for bivariate polynomials in (u, v) of total degree <= D
inline void d34e_enet_pmul(const double A[5][5], const double B[5][5], int D,
                           double C[5][5]) {
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j) C[i][j] = 0.0;
    for (int i = 0; i <= D; ++i)
        for (int j = 0; j <= D - i; ++j) {
            if (A[i][j] == 0.0) continue;
            for (int k = 0; k <= D - i; ++k)
                for (int l = 0; l <= D - i - j - k; ++l)
                    C[i + k][j + l] += A[i][j] * B[k][l];
        }
}

// x < 15: Taylor composition. o[k] is the derivative of order D with D - k
// derivatives in lambda and k in alpha.
inline void d34e_enet_emm_taylor(double lam, double al, int D, double* o) {
    const double b = 1.0 - al;
    const double x0 = lam * al / std::sqrt(lam * b);
    double m[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    d34e_enet_moments(x0, 6, m);
    double d[5];
    d34e_enet_m2der(m, d);
    // x(lambda + u, alpha + v) = sqrt(lambda + u) (alpha + v) (b - v)^(-1/2)
    double L[5], H[5], G[5];
    double cf = 1.0, lp = std::sqrt(lam);
    for (int i = 0; i <= D; ++i) {
        if (i > 0) {
            cf *= (0.5 - (i - 1)) / i;
            lp /= lam;
        }
        L[i] = cf * lp;
    }
    H[0] = 1.0 / std::sqrt(b);
    for (int j = 1; j <= D; ++j) H[j] = H[j - 1] * ((0.5 + (j - 1)) / j) / b;
    for (int j = 0; j <= D; ++j) G[j] = al * H[j] + (j >= 1 ? H[j - 1] : 0.0);
    double dl[5][5], pw[5][5], tmp[5][5], ps[5][5], cp[5][5], pe[5][5];
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j) {
            dl[i][j] = 0.0;
            pw[i][j] = 0.0;
            ps[i][j] = 0.0;
            cp[i][j] = 0.0;
        }
    for (int i = 0; i <= D; ++i)
        for (int j = 0; j <= D - i; ++j)
            if (i + j >= 1) dl[i][j] = L[i] * G[j];
    pw[0][0] = 1.0;
    double kfact = 1.0;
    for (int k = 0; k <= D; ++k) {
        if (k > 0) kfact *= k;
        const double f = d[k] / kfact;
        for (int i = 0; i <= D; ++i)
            for (int j = 0; j <= D - i; ++j) ps[i][j] += f * pw[i][j];
        d34e_enet_pmul(pw, dl, D, tmp);
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j) pw[i][j] = tmp[i][j];
    }
    // c = (lambda + u)(b - v)
    cp[0][0] = lam * b;
    cp[1][0] = b;
    cp[0][1] = -lam;
    cp[1][1] = -1.0;
    d34e_enet_pmul(cp, ps, D, pe);
    for (int k = 0; k <= D; ++k) {
        const int i = D - k, j = k;
        // -lambda^2 alpha^2 - 2 lambda (1 - alpha), differentiated exactly
        double poly = -d34e_enet_ff(2.0, i) * (i <= 2 ? std::pow(lam, 2 - i) : 0.0) *
                      d34e_enet_ff(2.0, j) * (j <= 2 ? std::pow(al, 2 - j) : 0.0);
        if (j == 0 && i <= 1)
            poly -= 2.0 * d34e_enet_ff(1.0, i) * std::pow(lam, 1 - i) * b;
        else if (j == 1 && i <= 1)
            poly += 2.0 * d34e_enet_ff(1.0, i) * std::pow(lam, 1 - i);
        double fi = 1.0, fj = 1.0;
        for (int q = 2; q <= i; ++q) fi *= q;
        for (int q = 2; q <= j; ++q) fj *= q;
        o[k] = poly + pe[i][j] * fi * fj;
    }
}

// x >= 15: the series, each monomial lambda^(2-n) b^n alpha^(2-2n) = a^2 t^n
// differentiated exactly (b = 1 - alpha)
inline void d34e_enet_emm_series(double lam, double al, int D, double* o) {
    static const double phi[31] = {
        1, 2, -2, 10, -74, 706, -8162, 110410, -1708394, 29752066, -576037442, 12277827850, -285764591114, 7213364729026, -196316804255522, 5731249477826890, -1.7867678947312182e+17, 5.9250857445438372e+18, -2.0825680275889236e+20, 7.7341580859426787e+21, -3.0262506729515713e+23, 1.2444051435099602e+25, -5.3650434562765592e+26, 2.4200145742440098e+28, -1.139863601904254e+30, 5.5963334851995967e+31, -2.8592787312798988e+33, 1.5179308863692238e+35, -8.3613927358484687e+36, 4.7727423221227011e+38, -2.8196190697395067e+40};
    const int NS = 30;
    const double b = 1.0 - al;
    const double t = b / (lam * al * al);
    const double a2 = lam * lam * al * al;
    for (int k = 0; k <= D; ++k) o[k] = 0.0;
    double tn = 1.0;
    for (int n = 0; n <= NS; ++n) {
        const double mn = a2 * tn;
        const double p = 2.0 - n, q = 2.0 - 2.0 * n;
        for (int k = 0; k <= D; ++k) {
            const int i = D - k, j = k;
            const double fi = d34e_enet_ff(p, i) / std::pow(lam, i);
            double s = 0.0, binom = 1.0;  // binom = choose(j, r)
            for (int r = 0; r <= j && r <= n; ++r) {
                if (r > 0) binom = binom * (j - r + 1) / r;
                s += binom * ((r % 2) ? -1.0 : 1.0) * d34e_enet_ff(n, r) /
                     std::pow(b, r) * d34e_enet_ff(q, j - r) /
                     std::pow(al, j - r);
            }
            o[k] -= phi[n] * mn * fi * s;
        }
        tn *= t;
    }
}

// the derivatives of order D of E_mm for the pair mu_mu: o[0..D]
inline void d34e_enet_emm(double lam, double al, int D, double* o) {
    const double x = lam * al / std::sqrt(lam * (1.0 - al));
    if (x < 15.0)
        d34e_enet_emm_taylor(lam, al, D, o);
    else
        d34e_enet_emm_series(lam, al, D, o);
}

}  // namespace


// The derivatives of order 3 of the expected information at one distinct
// pair (lambda, alpha): the 60 components in the order of d3expected_names().
inline void d34e_enet_comp3(double lam_u, double al_u, double U_u,
                                    const double* tab_u, double* o) {
    const double lam = lam_u;
    const double al = al_u;
    (void)U_u;  // the pair mu_mu does not read the table
    const double* Dt = tab_u;
    const double D03 = Dt[0 * 7 + 3];
    const double D04 = Dt[0 * 7 + 4];
    const double D05 = Dt[0 * 7 + 5];
    const double D12 = Dt[1 * 7 + 2];
    const double D13 = Dt[1 * 7 + 3];
    const double D14 = Dt[1 * 7 + 4];
    const double D21 = Dt[2 * 7 + 1];
    const double D22 = Dt[2 * 7 + 2];
    const double D23 = Dt[2 * 7 + 3];
    const double D30 = Dt[3 * 7 + 0];
    const double D31 = Dt[3 * 7 + 1];
    const double D32 = Dt[3 * 7 + 2];
    const double D40 = Dt[4 * 7 + 0];
    const double D41 = Dt[4 * 7 + 1];
    const double D50 = Dt[5 * 7 + 0];
        // mu_mu_mu_mu_mu
        o[0] = 0.0;
        // mu_mu_mu_mu_lambda
        o[1] = 0.0;
        // mu_mu_mu_mu_alpha
        o[2] = 0.0;
        // mu_mu_mu_lambda_lambda
        o[3] = 0.0;
        // mu_mu_mu_lambda_alpha
        o[4] = 0.0;
        // mu_mu_mu_alpha_alpha
        o[5] = 0.0;
        // mu_mu_(lambda, alpha) derivatives: E_mm = -c - 2 a U, computed apart from the table
        d34e_enet_emm(lam, al, 3, o + 6);
        // lambda_lambda_mu_mu_mu
        o[10] = 0.0;
        // lambda_lambda_mu_mu_lambda
        o[11] = 0.0;
        // lambda_lambda_mu_mu_alpha
        o[12] = 0.0;
        // lambda_lambda_mu_lambda_lambda
        o[13] = 0.0;
        // lambda_lambda_mu_lambda_alpha
        o[14] = 0.0;
        // lambda_lambda_mu_alpha_alpha
        o[15] = 0.0;
        // lambda_lambda_lambda_lambda_lambda
        o[16] = (((al - 1)*(al - 1)*(al - 1)*(al - 1)*(al - 1)))*(D05) + (-((al)*(al)*(al)*(al)*(al)))*(D50) + (-5*al*((al - 1)*(al - 1)*(al - 1)*(al - 1)))*(D14) + (10*((al)*(al))*((al - 1)*(al - 1)*(al - 1)))*(D23) + (-10*((al)*(al)*(al))*((al - 1)*(al - 1)))*(D32) + (5*((al)*(al)*(al)*(al))*(al - 1))*(D41);
        // lambda_lambda_lambda_lambda_alpha
        o[17] = (4*((al - 1)*(al - 1))*(4*al - 1))*(D13) + (-4*((al - 1)*(al - 1)*(al - 1)))*(D04) + (lam*((al - 1)*(al - 1)*(al - 1)*(al - 1)))*(D05) + (-lam*((al - 1)*(al - 1)*(al - 1))*(5*al - 1))*(D14) + (-12*al*(al - 1)*(2*al - 1))*(D22) + (4*((al)*(al))*(4*al - 3))*(D31) + (-4*((al)*(al)*(al)))*(D40) + (-((al)*(al)*(al)*(al))*lam)*(D50) + (2*al*lam*((al - 1)*(al - 1))*(5*al - 2))*(D23) + (-2*((al)*(al))*lam*(al - 1)*(5*al - 3))*(D32) + (((al)*(al)*(al))*lam*(5*al - 4))*(D41);
        // lambda_lambda_lambda_alpha_alpha
        o[18] = (6*(al - 1))*(D03) + (6*(3*al - 1))*(D21) + (-6*(3*al - 2))*(D12) + (((lam)*(lam))*((al - 1)*(al - 1)*(al - 1)))*(D05) + (((lam)*(lam))*(al - 1)*(10*((al)*(al)) - 8*al + 1))*(D23) + (-6*lam*((al - 1)*(al - 1)))*(D04) + (-6*lam*(6*((al)*(al)) - 6*al + 1))*(D22) + (-6*al)*(D30) + (-((lam)*(lam))*((al - 1)*(al - 1))*(5*al - 2))*(D14) + (12*lam*(al - 1)*(2*al - 1))*(D13) + (-((al)*(al)*(al))*((lam)*(lam)))*(D50) + (12*al*lam*(2*al - 1))*(D31) + (-al*((lam)*(lam))*(10*((al)*(al)) - 12*al + 3))*(D32) + (-6*((al)*(al))*lam)*(D40) + (((al)*(al))*((lam)*(lam))*(5*al - 3))*(D41);
        // lambda_lambda_alpha_alpha_alpha
        o[19] = (((lam)*(lam)*(lam))*((al - 1)*(al - 1)))*(D05) + (-((lam)*(lam)*(lam))*(10*((al)*(al)) - 8*al + 1))*(D32) + (-18*lam)*(D12) + (6*((lam)*(lam))*(4*al - 3))*(D13) + (-6*lam)*(D30) + (6*((lam)*(lam))*(4*al - 1))*(D31) + (-((lam)*(lam)*(lam))*(al - 1)*(5*al - 3))*(D14) + (((lam)*(lam)*(lam))*(10*((al)*(al)) - 12*al + 3))*(D23) + (6*lam)*(D03) + (-6*((lam)*(lam))*(al - 1))*(D04) + (18*lam)*(D21) + (-18*((lam)*(lam))*(2*al - 1))*(D22) + (-((al)*(al))*((lam)*(lam)*(lam)))*(D50) + (-6*al*((lam)*(lam)))*(D40) + (al*((lam)*(lam)*(lam))*(5*al - 2))*(D41);
        // alpha_alpha_mu_mu_mu
        o[20] = 0.0;
        // alpha_alpha_mu_mu_lambda
        o[21] = 0.0;
        // alpha_alpha_mu_mu_alpha
        o[22] = 0.0;
        // alpha_alpha_mu_lambda_lambda
        o[23] = 0.0;
        // alpha_alpha_mu_lambda_alpha
        o[24] = 0.0;
        // alpha_alpha_mu_alpha_alpha
        o[25] = 0.0;
        // alpha_alpha_lambda_lambda_lambda
        o[26] = (6*(al - 1))*(D03) + (6*(3*al - 1))*(D21) + (-6*(3*al - 2))*(D12) + (((lam)*(lam))*((al - 1)*(al - 1)*(al - 1)))*(D05) + (((lam)*(lam))*(al - 1)*(10*((al)*(al)) - 8*al + 1))*(D23) + (-6*lam*((al - 1)*(al - 1)))*(D04) + (-6*lam*(6*((al)*(al)) - 6*al + 1))*(D22) + (-6*al)*(D30) + (-((lam)*(lam))*((al - 1)*(al - 1))*(5*al - 2))*(D14) + (12*lam*(al - 1)*(2*al - 1))*(D13) + (-((al)*(al)*(al))*((lam)*(lam)))*(D50) + (12*al*lam*(2*al - 1))*(D31) + (-al*((lam)*(lam))*(10*((al)*(al)) - 12*al + 3))*(D32) + (-6*((al)*(al))*lam)*(D40) + (((al)*(al))*((lam)*(lam))*(5*al - 3))*(D41);
        // alpha_alpha_lambda_lambda_alpha
        o[27] = (((lam)*(lam)*(lam))*((al - 1)*(al - 1)))*(D05) + (-((lam)*(lam)*(lam))*(10*((al)*(al)) - 8*al + 1))*(D32) + (-18*lam)*(D12) + (6*((lam)*(lam))*(4*al - 3))*(D13) + (-6*lam)*(D30) + (6*((lam)*(lam))*(4*al - 1))*(D31) + (-((lam)*(lam)*(lam))*(al - 1)*(5*al - 3))*(D14) + (((lam)*(lam)*(lam))*(10*((al)*(al)) - 12*al + 3))*(D23) + (6*lam)*(D03) + (-6*((lam)*(lam))*(al - 1))*(D04) + (18*lam)*(D21) + (-18*((lam)*(lam))*(2*al - 1))*(D22) + (-((al)*(al))*((lam)*(lam)*(lam)))*(D50) + (-6*al*((lam)*(lam)))*(D40) + (al*((lam)*(lam)*(lam))*(5*al - 2))*(D41);
        // alpha_alpha_lambda_alpha_alpha
        o[28] = (((lam)*(lam)*(lam)*(lam))*(al - 1))*(D05) + (((lam)*(lam)*(lam)*(lam))*(5*al - 1))*(D41) + (-24*((lam)*(lam)*(lam)))*(D22) + (2*((lam)*(lam)*(lam)*(lam))*(5*al - 3))*(D23) + (-4*((lam)*(lam)*(lam)))*(D04) + (-4*((lam)*(lam)*(lam)))*(D40) + (-((lam)*(lam)*(lam)*(lam))*(5*al - 4))*(D14) + (-2*((lam)*(lam)*(lam)*(lam))*(5*al - 2))*(D32) + (16*((lam)*(lam)*(lam)))*(D13) + (16*((lam)*(lam)*(lam)))*(D31) + (-al*((lam)*(lam)*(lam)*(lam)))*(D50);
        // alpha_alpha_alpha_alpha_alpha
        o[29] = (((lam)*(lam)*(lam)*(lam)*(lam)))*(D05) + (-((lam)*(lam)*(lam)*(lam)*(lam)))*(D50) + (-10*((lam)*(lam)*(lam)*(lam)*(lam)))*(D32) + (-5*((lam)*(lam)*(lam)*(lam)*(lam)))*(D14) + (5*((lam)*(lam)*(lam)*(lam)*(lam)))*(D41) + (10*((lam)*(lam)*(lam)*(lam)*(lam)))*(D23);
        // mu_lambda_mu_mu_mu
        o[30] = 0.0;
        // mu_lambda_mu_mu_lambda
        o[31] = 0.0;
        // mu_lambda_mu_mu_alpha
        o[32] = 0.0;
        // mu_lambda_mu_lambda_lambda
        o[33] = 0.0;
        // mu_lambda_mu_lambda_alpha
        o[34] = 0.0;
        // mu_lambda_mu_alpha_alpha
        o[35] = 0.0;
        // mu_lambda_lambda_lambda_lambda
        o[36] = 0.0;
        // mu_lambda_lambda_lambda_alpha
        o[37] = 0.0;
        // mu_lambda_lambda_alpha_alpha
        o[38] = 0.0;
        // mu_lambda_alpha_alpha_alpha
        o[39] = 0.0;
        // mu_alpha_mu_mu_mu
        o[40] = 0.0;
        // mu_alpha_mu_mu_lambda
        o[41] = 0.0;
        // mu_alpha_mu_mu_alpha
        o[42] = 0.0;
        // mu_alpha_mu_lambda_lambda
        o[43] = 0.0;
        // mu_alpha_mu_lambda_alpha
        o[44] = 0.0;
        // mu_alpha_mu_alpha_alpha
        o[45] = 0.0;
        // mu_alpha_lambda_lambda_lambda
        o[46] = 0.0;
        // mu_alpha_lambda_lambda_alpha
        o[47] = 0.0;
        // mu_alpha_lambda_alpha_alpha
        o[48] = 0.0;
        // mu_alpha_alpha_alpha_alpha
        o[49] = 0.0;
        // lambda_alpha_mu_mu_mu
        o[50] = 0.0;
        // lambda_alpha_mu_mu_lambda
        o[51] = 0.0;
        // lambda_alpha_mu_mu_alpha
        o[52] = 0.0;
        // lambda_alpha_mu_lambda_lambda
        o[53] = 0.0;
        // lambda_alpha_mu_lambda_alpha
        o[54] = 0.0;
        // lambda_alpha_mu_alpha_alpha
        o[55] = 0.0;
        // lambda_alpha_lambda_lambda_lambda
        o[56] = (3*((al - 1)*(al - 1))*(4*al - 1))*(D13) + (-3*((al - 1)*(al - 1)*(al - 1)))*(D04) + (lam*((al - 1)*(al - 1)*(al - 1)*(al - 1)))*(D05) + (-lam*((al - 1)*(al - 1)*(al - 1))*(5*al - 1))*(D14) + (-9*al*(al - 1)*(2*al - 1))*(D22) + (3*((al)*(al))*(4*al - 3))*(D31) + (-3*((al)*(al)*(al)))*(D40) + (-((al)*(al)*(al)*(al))*lam)*(D50) + (2*al*lam*((al - 1)*(al - 1))*(5*al - 2))*(D23) + (-2*((al)*(al))*lam*(al - 1)*(5*al - 3))*(D32) + (((al)*(al)*(al))*lam*(5*al - 4))*(D41);
        // lambda_alpha_lambda_lambda_alpha
        o[57] = (4*(al - 1))*(D03) + (4*(3*al - 1))*(D21) + (-4*(3*al - 2))*(D12) + (((lam)*(lam))*((al - 1)*(al - 1)*(al - 1)))*(D05) + (((lam)*(lam))*(al - 1)*(10*((al)*(al)) - 8*al + 1))*(D23) + (-5*lam*((al - 1)*(al - 1)))*(D04) + (-5*lam*(6*((al)*(al)) - 6*al + 1))*(D22) + (-4*al)*(D30) + (-((lam)*(lam))*((al - 1)*(al - 1))*(5*al - 2))*(D14) + (10*lam*(al - 1)*(2*al - 1))*(D13) + (-((al)*(al)*(al))*((lam)*(lam)))*(D50) + (10*al*lam*(2*al - 1))*(D31) + (-al*((lam)*(lam))*(10*((al)*(al)) - 12*al + 3))*(D32) + (-5*((al)*(al))*lam)*(D40) + (((al)*(al))*((lam)*(lam))*(5*al - 3))*(D41);
        // lambda_alpha_lambda_alpha_alpha
        o[58] = (((lam)*(lam)*(lam))*((al - 1)*(al - 1)))*(D05) + (-((lam)*(lam)*(lam))*(10*((al)*(al)) - 8*al + 1))*(D32) + (5*((lam)*(lam))*(4*al - 3))*(D13) + (-12*lam)*(D12) + (5*((lam)*(lam))*(4*al - 1))*(D31) + (-4*lam)*(D30) + (-((lam)*(lam)*(lam))*(al - 1)*(5*al - 3))*(D14) + (((lam)*(lam)*(lam))*(10*((al)*(al)) - 12*al + 3))*(D23) + (4*lam)*(D03) + (-5*((lam)*(lam))*(al - 1))*(D04) + (12*lam)*(D21) + (-15*((lam)*(lam))*(2*al - 1))*(D22) + (-((al)*(al))*((lam)*(lam)*(lam)))*(D50) + (-5*al*((lam)*(lam)))*(D40) + (al*((lam)*(lam)*(lam))*(5*al - 2))*(D41);
        // lambda_alpha_alpha_alpha_alpha
        o[59] = (((lam)*(lam)*(lam)*(lam))*(al - 1))*(D05) + (((lam)*(lam)*(lam)*(lam))*(5*al - 1))*(D41) + (-18*((lam)*(lam)*(lam)))*(D22) + (2*((lam)*(lam)*(lam)*(lam))*(5*al - 3))*(D23) + (-3*((lam)*(lam)*(lam)))*(D04) + (-3*((lam)*(lam)*(lam)))*(D40) + (-((lam)*(lam)*(lam)*(lam))*(5*al - 4))*(D14) + (-2*((lam)*(lam)*(lam)*(lam))*(5*al - 2))*(D32) + (12*((lam)*(lam)*(lam)))*(D13) + (12*((lam)*(lam)*(lam)))*(D31) + (-al*((lam)*(lam)*(lam)*(lam)))*(D50);
}


// The derivatives of order 4 of the expected information at one distinct
// pair (lambda, alpha): the 90 components in the order of d4expected_names().
inline void d34e_enet_comp4(double lam_u, double al_u, double U_u,
                                    const double* tab_u, double* o) {
    const double lam = lam_u;
    const double al = al_u;
    (void)U_u;  // the pair mu_mu does not read the table
    const double* Dt = tab_u;
    const double D03 = Dt[0 * 7 + 3];
    const double D04 = Dt[0 * 7 + 4];
    const double D05 = Dt[0 * 7 + 5];
    const double D06 = Dt[0 * 7 + 6];
    const double D12 = Dt[1 * 7 + 2];
    const double D13 = Dt[1 * 7 + 3];
    const double D14 = Dt[1 * 7 + 4];
    const double D15 = Dt[1 * 7 + 5];
    const double D21 = Dt[2 * 7 + 1];
    const double D22 = Dt[2 * 7 + 2];
    const double D23 = Dt[2 * 7 + 3];
    const double D24 = Dt[2 * 7 + 4];
    const double D30 = Dt[3 * 7 + 0];
    const double D31 = Dt[3 * 7 + 1];
    const double D32 = Dt[3 * 7 + 2];
    const double D33 = Dt[3 * 7 + 3];
    const double D40 = Dt[4 * 7 + 0];
    const double D41 = Dt[4 * 7 + 1];
    const double D42 = Dt[4 * 7 + 2];
    const double D50 = Dt[5 * 7 + 0];
    const double D51 = Dt[5 * 7 + 1];
    const double D60 = Dt[6 * 7 + 0];
        // mu_mu_mu_mu_mu_mu
        o[0] = 0.0;
        // mu_mu_mu_mu_mu_lambda
        o[1] = 0.0;
        // mu_mu_mu_mu_mu_alpha
        o[2] = 0.0;
        // mu_mu_mu_mu_lambda_lambda
        o[3] = 0.0;
        // mu_mu_mu_mu_lambda_alpha
        o[4] = 0.0;
        // mu_mu_mu_mu_alpha_alpha
        o[5] = 0.0;
        // mu_mu_mu_lambda_lambda_lambda
        o[6] = 0.0;
        // mu_mu_mu_lambda_lambda_alpha
        o[7] = 0.0;
        // mu_mu_mu_lambda_alpha_alpha
        o[8] = 0.0;
        // mu_mu_mu_alpha_alpha_alpha
        o[9] = 0.0;
        // mu_mu_(lambda, alpha) derivatives: E_mm = -c - 2 a U, computed apart from the table
        d34e_enet_emm(lam, al, 4, o + 10);
        // lambda_lambda_mu_mu_mu_mu
        o[15] = 0.0;
        // lambda_lambda_mu_mu_mu_lambda
        o[16] = 0.0;
        // lambda_lambda_mu_mu_mu_alpha
        o[17] = 0.0;
        // lambda_lambda_mu_mu_lambda_lambda
        o[18] = 0.0;
        // lambda_lambda_mu_mu_lambda_alpha
        o[19] = 0.0;
        // lambda_lambda_mu_mu_alpha_alpha
        o[20] = 0.0;
        // lambda_lambda_mu_lambda_lambda_lambda
        o[21] = 0.0;
        // lambda_lambda_mu_lambda_lambda_alpha
        o[22] = 0.0;
        // lambda_lambda_mu_lambda_alpha_alpha
        o[23] = 0.0;
        // lambda_lambda_mu_alpha_alpha_alpha
        o[24] = 0.0;
        // lambda_lambda_lambda_lambda_lambda_lambda
        o[25] = (-((al - 1)*(al - 1)*(al - 1)*(al - 1)*(al - 1)*(al - 1)))*(D06) + (-((al)*(al)*(al)*(al)*(al)*(al)))*(D60) + (-15*((al)*(al))*((al - 1)*(al - 1)*(al - 1)*(al - 1)))*(D24) + (6*al*((al - 1)*(al - 1)*(al - 1)*(al - 1)*(al - 1)))*(D15) + (20*((al)*(al)*(al))*((al - 1)*(al - 1)*(al - 1)))*(D33) + (-15*((al)*(al)*(al)*(al))*((al - 1)*(al - 1)))*(D42) + (6*((al)*(al)*(al)*(al)*(al))*(al - 1))*(D51);
        // lambda_lambda_lambda_lambda_lambda_alpha
        o[26] = (-5*((al - 1)*(al - 1)*(al - 1))*(5*al - 1))*(D14) + (5*((al - 1)*(al - 1)*(al - 1)*(al - 1)))*(D05) + (-lam*((al - 1)*(al - 1)*(al - 1)*(al - 1)*(al - 1)))*(D06) + (lam*((al - 1)*(al - 1)*(al - 1)*(al - 1))*(6*al - 1))*(D15) + (10*al*((al - 1)*(al - 1))*(5*al - 2))*(D23) + (-10*((al)*(al))*(al - 1)*(5*al - 3))*(D32) + (5*((al)*(al)*(al))*(5*al - 4))*(D41) + (-5*((al)*(al)*(al)*(al)))*(D50) + (-((al)*(al)*(al)*(al)*(al))*lam)*(D60) + (-5*al*lam*((al - 1)*(al - 1)*(al - 1))*(3*al - 1))*(D24) + (10*((al)*(al))*lam*((al - 1)*(al - 1))*(2*al - 1))*(D33) + (-5*((al)*(al)*(al))*lam*(al - 1)*(3*al - 2))*(D42) + (((al)*(al)*(al)*(al))*lam*(6*al - 5))*(D51);
        // lambda_lambda_lambda_lambda_alpha_alpha
        o[27] = (-12*((al - 1)*(al - 1)))*(D04) + (-12*(6*((al)*(al)) - 6*al + 1))*(D22) + (24*(al - 1)*(2*al - 1))*(D13) + (-((lam)*(lam))*((al - 1)*(al - 1)*(al - 1)*(al - 1)))*(D06) + (-((lam)*(lam))*((al - 1)*(al - 1))*(15*((al)*(al)) - 10*al + 1))*(D24) + (24*al*(2*al - 1))*(D31) + (-12*((al)*(al)))*(D40) + (8*lam*((al - 1)*(al - 1)*(al - 1)))*(D05) + (8*lam*(al - 1)*(10*((al)*(al)) - 8*al + 1))*(D23) + (2*((lam)*(lam))*((al - 1)*(al - 1)*(al - 1))*(3*al - 1))*(D15) + (-8*lam*((al - 1)*(al - 1))*(5*al - 2))*(D14) + (-((al)*(al)*(al)*(al))*((lam)*(lam)))*(D60) + (-8*al*lam*(10*((al)*(al)) - 12*al + 3))*(D32) + (4*al*((lam)*(lam))*(al - 1)*(5*((al)*(al)) - 5*al + 1))*(D33) + (8*((al)*(al))*lam*(5*al - 3))*(D41) + (-((al)*(al))*((lam)*(lam))*(15*((al)*(al)) - 20*al + 6))*(D42) + (-8*((al)*(al)*(al))*lam)*(D50) + (2*((al)*(al)*(al))*((lam)*(lam))*(3*al - 2))*(D51);
        // lambda_lambda_lambda_alpha_alpha_alpha
        o[28] = (-18)*(D12) + (-6)*(D30) + (6)*(D03) + (18)*(D21) + (-((lam)*(lam)*(lam))*((al - 1)*(al - 1)*(al - 1)))*(D06) + (((lam)*(lam)*(lam))*(2*al - 1)*(10*((al)*(al)) - 10*al + 1))*(D33) + (18*lam*(4*al - 3))*(D13) + (-9*((lam)*(lam))*(al - 1)*(5*al - 3))*(D14) + (18*lam*(4*al - 1))*(D31) + (-9*((lam)*(lam))*(10*((al)*(al)) - 8*al + 1))*(D32) + (3*((lam)*(lam)*(lam))*((al - 1)*(al - 1))*(2*al - 1))*(D15) + (-3*((lam)*(lam)*(lam))*(al - 1)*(5*((al)*(al)) - 5*al + 1))*(D24) + (9*((lam)*(lam))*((al - 1)*(al - 1)))*(D05) + (-18*lam*(al - 1))*(D04) + (9*((lam)*(lam))*(10*((al)*(al)) - 12*al + 3))*(D23) + (-54*lam*(2*al - 1))*(D22) + (-((al)*(al)*(al))*((lam)*(lam)*(lam)))*(D60) + (-18*al*lam)*(D40) + (9*al*((lam)*(lam))*(5*al - 2))*(D41) + (-3*al*((lam)*(lam)*(lam))*(5*((al)*(al)) - 5*al + 1))*(D42) + (-9*((al)*(al))*((lam)*(lam)))*(D50) + (3*((al)*(al))*((lam)*(lam)*(lam))*(2*al - 1))*(D51);
        // lambda_lambda_alpha_alpha_alpha_alpha
        o[29] = (-((lam)*(lam)*(lam)*(lam))*((al - 1)*(al - 1)))*(D06) + (-((lam)*(lam)*(lam)*(lam))*(15*((al)*(al)) - 10*al + 1))*(D42) + (-72*((lam)*(lam)))*(D22) + (16*((lam)*(lam)*(lam))*(5*al - 3))*(D23) + (-12*((lam)*(lam)))*(D04) + (-12*((lam)*(lam)))*(D40) + (8*((lam)*(lam)*(lam))*(al - 1))*(D05) + (8*((lam)*(lam)*(lam))*(5*al - 1))*(D41) + (-((lam)*(lam)*(lam)*(lam))*(15*((al)*(al)) - 20*al + 6))*(D24) + (2*((lam)*(lam)*(lam)*(lam))*(al - 1)*(3*al - 2))*(D15) + (4*((lam)*(lam)*(lam)*(lam))*(5*((al)*(al)) - 5*al + 1))*(D33) + (-8*((lam)*(lam)*(lam))*(5*al - 4))*(D14) + (-16*((lam)*(lam)*(lam))*(5*al - 2))*(D32) + (48*((lam)*(lam)))*(D13) + (48*((lam)*(lam)))*(D31) + (-((al)*(al))*((lam)*(lam)*(lam)*(lam)))*(D60) + (-8*al*((lam)*(lam)*(lam)))*(D50) + (2*al*((lam)*(lam)*(lam)*(lam))*(3*al - 1))*(D51);
        // alpha_alpha_mu_mu_mu_mu
        o[30] = 0.0;
        // alpha_alpha_mu_mu_mu_lambda
        o[31] = 0.0;
        // alpha_alpha_mu_mu_mu_alpha
        o[32] = 0.0;
        // alpha_alpha_mu_mu_lambda_lambda
        o[33] = 0.0;
        // alpha_alpha_mu_mu_lambda_alpha
        o[34] = 0.0;
        // alpha_alpha_mu_mu_alpha_alpha
        o[35] = 0.0;
        // alpha_alpha_mu_lambda_lambda_lambda
        o[36] = 0.0;
        // alpha_alpha_mu_lambda_lambda_alpha
        o[37] = 0.0;
        // alpha_alpha_mu_lambda_alpha_alpha
        o[38] = 0.0;
        // alpha_alpha_mu_alpha_alpha_alpha
        o[39] = 0.0;
        // alpha_alpha_lambda_lambda_lambda_lambda
        o[40] = (-12*((al - 1)*(al - 1)))*(D04) + (-12*(6*((al)*(al)) - 6*al + 1))*(D22) + (24*(al - 1)*(2*al - 1))*(D13) + (-((lam)*(lam))*((al - 1)*(al - 1)*(al - 1)*(al - 1)))*(D06) + (-((lam)*(lam))*((al - 1)*(al - 1))*(15*((al)*(al)) - 10*al + 1))*(D24) + (24*al*(2*al - 1))*(D31) + (-12*((al)*(al)))*(D40) + (8*lam*((al - 1)*(al - 1)*(al - 1)))*(D05) + (8*lam*(al - 1)*(10*((al)*(al)) - 8*al + 1))*(D23) + (2*((lam)*(lam))*((al - 1)*(al - 1)*(al - 1))*(3*al - 1))*(D15) + (-8*lam*((al - 1)*(al - 1))*(5*al - 2))*(D14) + (-((al)*(al)*(al)*(al))*((lam)*(lam)))*(D60) + (-8*al*lam*(10*((al)*(al)) - 12*al + 3))*(D32) + (4*al*((lam)*(lam))*(al - 1)*(5*((al)*(al)) - 5*al + 1))*(D33) + (8*((al)*(al))*lam*(5*al - 3))*(D41) + (-((al)*(al))*((lam)*(lam))*(15*((al)*(al)) - 20*al + 6))*(D42) + (-8*((al)*(al)*(al))*lam)*(D50) + (2*((al)*(al)*(al))*((lam)*(lam))*(3*al - 2))*(D51);
        // alpha_alpha_lambda_lambda_lambda_alpha
        o[41] = (-18)*(D12) + (-6)*(D30) + (6)*(D03) + (18)*(D21) + (-((lam)*(lam)*(lam))*((al - 1)*(al - 1)*(al - 1)))*(D06) + (((lam)*(lam)*(lam))*(2*al - 1)*(10*((al)*(al)) - 10*al + 1))*(D33) + (18*lam*(4*al - 3))*(D13) + (-9*((lam)*(lam))*(al - 1)*(5*al - 3))*(D14) + (18*lam*(4*al - 1))*(D31) + (-9*((lam)*(lam))*(10*((al)*(al)) - 8*al + 1))*(D32) + (3*((lam)*(lam)*(lam))*((al - 1)*(al - 1))*(2*al - 1))*(D15) + (-3*((lam)*(lam)*(lam))*(al - 1)*(5*((al)*(al)) - 5*al + 1))*(D24) + (9*((lam)*(lam))*((al - 1)*(al - 1)))*(D05) + (-18*lam*(al - 1))*(D04) + (9*((lam)*(lam))*(10*((al)*(al)) - 12*al + 3))*(D23) + (-54*lam*(2*al - 1))*(D22) + (-((al)*(al)*(al))*((lam)*(lam)*(lam)))*(D60) + (-18*al*lam)*(D40) + (9*al*((lam)*(lam))*(5*al - 2))*(D41) + (-3*al*((lam)*(lam)*(lam))*(5*((al)*(al)) - 5*al + 1))*(D42) + (-9*((al)*(al))*((lam)*(lam)))*(D50) + (3*((al)*(al))*((lam)*(lam)*(lam))*(2*al - 1))*(D51);
        // alpha_alpha_lambda_lambda_alpha_alpha
        o[42] = (-((lam)*(lam)*(lam)*(lam))*((al - 1)*(al - 1)))*(D06) + (-((lam)*(lam)*(lam)*(lam))*(15*((al)*(al)) - 10*al + 1))*(D42) + (-72*((lam)*(lam)))*(D22) + (16*((lam)*(lam)*(lam))*(5*al - 3))*(D23) + (-12*((lam)*(lam)))*(D04) + (-12*((lam)*(lam)))*(D40) + (8*((lam)*(lam)*(lam))*(al - 1))*(D05) + (8*((lam)*(lam)*(lam))*(5*al - 1))*(D41) + (-((lam)*(lam)*(lam)*(lam))*(15*((al)*(al)) - 20*al + 6))*(D24) + (2*((lam)*(lam)*(lam)*(lam))*(al - 1)*(3*al - 2))*(D15) + (4*((lam)*(lam)*(lam)*(lam))*(5*((al)*(al)) - 5*al + 1))*(D33) + (-8*((lam)*(lam)*(lam))*(5*al - 4))*(D14) + (-16*((lam)*(lam)*(lam))*(5*al - 2))*(D32) + (48*((lam)*(lam)))*(D13) + (48*((lam)*(lam)))*(D31) + (-((al)*(al))*((lam)*(lam)*(lam)*(lam)))*(D60) + (-8*al*((lam)*(lam)*(lam)))*(D50) + (2*al*((lam)*(lam)*(lam)*(lam))*(3*al - 1))*(D51);
        // alpha_alpha_lambda_alpha_alpha_alpha
        o[43] = (-((lam)*(lam)*(lam)*(lam)*(lam))*(al - 1))*(D06) + (((lam)*(lam)*(lam)*(lam)*(lam))*(6*al - 1))*(D51) + (-50*((lam)*(lam)*(lam)*(lam)))*(D32) + (-25*((lam)*(lam)*(lam)*(lam)))*(D14) + (10*((lam)*(lam)*(lam)*(lam)*(lam))*(2*al - 1))*(D33) + (((lam)*(lam)*(lam)*(lam)*(lam))*(6*al - 5))*(D15) + (-5*((lam)*(lam)*(lam)*(lam)))*(D50) + (5*((lam)*(lam)*(lam)*(lam)))*(D05) + (-5*((lam)*(lam)*(lam)*(lam)*(lam))*(3*al - 1))*(D42) + (-5*((lam)*(lam)*(lam)*(lam)*(lam))*(3*al - 2))*(D24) + (25*((lam)*(lam)*(lam)*(lam)))*(D41) + (50*((lam)*(lam)*(lam)*(lam)))*(D23) + (-al*((lam)*(lam)*(lam)*(lam)*(lam)))*(D60);
        // alpha_alpha_alpha_alpha_alpha_alpha
        o[44] = (-((lam)*(lam)*(lam)*(lam)*(lam)*(lam)))*(D06) + (-((lam)*(lam)*(lam)*(lam)*(lam)*(lam)))*(D60) + (-15*((lam)*(lam)*(lam)*(lam)*(lam)*(lam)))*(D24) + (-15*((lam)*(lam)*(lam)*(lam)*(lam)*(lam)))*(D42) + (6*((lam)*(lam)*(lam)*(lam)*(lam)*(lam)))*(D15) + (6*((lam)*(lam)*(lam)*(lam)*(lam)*(lam)))*(D51) + (20*((lam)*(lam)*(lam)*(lam)*(lam)*(lam)))*(D33);
        // mu_lambda_mu_mu_mu_mu
        o[45] = 0.0;
        // mu_lambda_mu_mu_mu_lambda
        o[46] = 0.0;
        // mu_lambda_mu_mu_mu_alpha
        o[47] = 0.0;
        // mu_lambda_mu_mu_lambda_lambda
        o[48] = 0.0;
        // mu_lambda_mu_mu_lambda_alpha
        o[49] = 0.0;
        // mu_lambda_mu_mu_alpha_alpha
        o[50] = 0.0;
        // mu_lambda_mu_lambda_lambda_lambda
        o[51] = 0.0;
        // mu_lambda_mu_lambda_lambda_alpha
        o[52] = 0.0;
        // mu_lambda_mu_lambda_alpha_alpha
        o[53] = 0.0;
        // mu_lambda_mu_alpha_alpha_alpha
        o[54] = 0.0;
        // mu_lambda_lambda_lambda_lambda_lambda
        o[55] = 0.0;
        // mu_lambda_lambda_lambda_lambda_alpha
        o[56] = 0.0;
        // mu_lambda_lambda_lambda_alpha_alpha
        o[57] = 0.0;
        // mu_lambda_lambda_alpha_alpha_alpha
        o[58] = 0.0;
        // mu_lambda_alpha_alpha_alpha_alpha
        o[59] = 0.0;
        // mu_alpha_mu_mu_mu_mu
        o[60] = 0.0;
        // mu_alpha_mu_mu_mu_lambda
        o[61] = 0.0;
        // mu_alpha_mu_mu_mu_alpha
        o[62] = 0.0;
        // mu_alpha_mu_mu_lambda_lambda
        o[63] = 0.0;
        // mu_alpha_mu_mu_lambda_alpha
        o[64] = 0.0;
        // mu_alpha_mu_mu_alpha_alpha
        o[65] = 0.0;
        // mu_alpha_mu_lambda_lambda_lambda
        o[66] = 0.0;
        // mu_alpha_mu_lambda_lambda_alpha
        o[67] = 0.0;
        // mu_alpha_mu_lambda_alpha_alpha
        o[68] = 0.0;
        // mu_alpha_mu_alpha_alpha_alpha
        o[69] = 0.0;
        // mu_alpha_lambda_lambda_lambda_lambda
        o[70] = 0.0;
        // mu_alpha_lambda_lambda_lambda_alpha
        o[71] = 0.0;
        // mu_alpha_lambda_lambda_alpha_alpha
        o[72] = 0.0;
        // mu_alpha_lambda_alpha_alpha_alpha
        o[73] = 0.0;
        // mu_alpha_alpha_alpha_alpha_alpha
        o[74] = 0.0;
        // lambda_alpha_mu_mu_mu_mu
        o[75] = 0.0;
        // lambda_alpha_mu_mu_mu_lambda
        o[76] = 0.0;
        // lambda_alpha_mu_mu_mu_alpha
        o[77] = 0.0;
        // lambda_alpha_mu_mu_lambda_lambda
        o[78] = 0.0;
        // lambda_alpha_mu_mu_lambda_alpha
        o[79] = 0.0;
        // lambda_alpha_mu_mu_alpha_alpha
        o[80] = 0.0;
        // lambda_alpha_mu_lambda_lambda_lambda
        o[81] = 0.0;
        // lambda_alpha_mu_lambda_lambda_alpha
        o[82] = 0.0;
        // lambda_alpha_mu_lambda_alpha_alpha
        o[83] = 0.0;
        // lambda_alpha_mu_alpha_alpha_alpha
        o[84] = 0.0;
        // lambda_alpha_lambda_lambda_lambda_lambda
        o[85] = (-4*((al - 1)*(al - 1)*(al - 1))*(5*al - 1))*(D14) + (4*((al - 1)*(al - 1)*(al - 1)*(al - 1)))*(D05) + (-lam*((al - 1)*(al - 1)*(al - 1)*(al - 1)*(al - 1)))*(D06) + (lam*((al - 1)*(al - 1)*(al - 1)*(al - 1))*(6*al - 1))*(D15) + (8*al*((al - 1)*(al - 1))*(5*al - 2))*(D23) + (-8*((al)*(al))*(al - 1)*(5*al - 3))*(D32) + (4*((al)*(al)*(al))*(5*al - 4))*(D41) + (-4*((al)*(al)*(al)*(al)))*(D50) + (-((al)*(al)*(al)*(al)*(al))*lam)*(D60) + (-5*al*lam*((al - 1)*(al - 1)*(al - 1))*(3*al - 1))*(D24) + (10*((al)*(al))*lam*((al - 1)*(al - 1))*(2*al - 1))*(D33) + (-5*((al)*(al)*(al))*lam*(al - 1)*(3*al - 2))*(D42) + (((al)*(al)*(al)*(al))*lam*(6*al - 5))*(D51);
        // lambda_alpha_lambda_lambda_lambda_alpha
        o[86] = (-9*((al - 1)*(al - 1)))*(D04) + (-9*(6*((al)*(al)) - 6*al + 1))*(D22) + (18*(al - 1)*(2*al - 1))*(D13) + (-((lam)*(lam))*((al - 1)*(al - 1)*(al - 1)*(al - 1)))*(D06) + (-((lam)*(lam))*((al - 1)*(al - 1))*(15*((al)*(al)) - 10*al + 1))*(D24) + (18*al*(2*al - 1))*(D31) + (-9*((al)*(al)))*(D40) + (7*lam*((al - 1)*(al - 1)*(al - 1)))*(D05) + (7*lam*(al - 1)*(10*((al)*(al)) - 8*al + 1))*(D23) + (2*((lam)*(lam))*((al - 1)*(al - 1)*(al - 1))*(3*al - 1))*(D15) + (-7*lam*((al - 1)*(al - 1))*(5*al - 2))*(D14) + (-((al)*(al)*(al)*(al))*((lam)*(lam)))*(D60) + (-7*al*lam*(10*((al)*(al)) - 12*al + 3))*(D32) + (4*al*((lam)*(lam))*(al - 1)*(5*((al)*(al)) - 5*al + 1))*(D33) + (7*((al)*(al))*lam*(5*al - 3))*(D41) + (-((al)*(al))*((lam)*(lam))*(15*((al)*(al)) - 20*al + 6))*(D42) + (-7*((al)*(al)*(al))*lam)*(D50) + (2*((al)*(al)*(al))*((lam)*(lam))*(3*al - 2))*(D51);
        // lambda_alpha_lambda_lambda_alpha_alpha
        o[87] = (-12)*(D12) + (-4)*(D30) + (4)*(D03) + (12)*(D21) + (-((lam)*(lam)*(lam))*((al - 1)*(al - 1)*(al - 1)))*(D06) + (((lam)*(lam)*(lam))*(2*al - 1)*(10*((al)*(al)) - 10*al + 1))*(D33) + (14*lam*(4*al - 3))*(D13) + (-8*((lam)*(lam))*(al - 1)*(5*al - 3))*(D14) + (14*lam*(4*al - 1))*(D31) + (-8*((lam)*(lam))*(10*((al)*(al)) - 8*al + 1))*(D32) + (3*((lam)*(lam)*(lam))*((al - 1)*(al - 1))*(2*al - 1))*(D15) + (-3*((lam)*(lam)*(lam))*(al - 1)*(5*((al)*(al)) - 5*al + 1))*(D24) + (8*((lam)*(lam))*((al - 1)*(al - 1)))*(D05) + (-14*lam*(al - 1))*(D04) + (8*((lam)*(lam))*(10*((al)*(al)) - 12*al + 3))*(D23) + (-42*lam*(2*al - 1))*(D22) + (-((al)*(al)*(al))*((lam)*(lam)*(lam)))*(D60) + (8*al*((lam)*(lam))*(5*al - 2))*(D41) + (-3*al*((lam)*(lam)*(lam))*(5*((al)*(al)) - 5*al + 1))*(D42) + (-14*al*lam)*(D40) + (-8*((al)*(al))*((lam)*(lam)))*(D50) + (3*((al)*(al))*((lam)*(lam)*(lam))*(2*al - 1))*(D51);
        // lambda_alpha_lambda_alpha_alpha_alpha
        o[88] = (-((lam)*(lam)*(lam)*(lam))*((al - 1)*(al - 1)))*(D06) + (-((lam)*(lam)*(lam)*(lam))*(15*((al)*(al)) - 10*al + 1))*(D42) + (-54*((lam)*(lam)))*(D22) + (14*((lam)*(lam)*(lam))*(5*al - 3))*(D23) + (-9*((lam)*(lam)))*(D04) + (-9*((lam)*(lam)))*(D40) + (7*((lam)*(lam)*(lam))*(al - 1))*(D05) + (7*((lam)*(lam)*(lam))*(5*al - 1))*(D41) + (-((lam)*(lam)*(lam)*(lam))*(15*((al)*(al)) - 20*al + 6))*(D24) + (2*((lam)*(lam)*(lam)*(lam))*(al - 1)*(3*al - 2))*(D15) + (4*((lam)*(lam)*(lam)*(lam))*(5*((al)*(al)) - 5*al + 1))*(D33) + (-7*((lam)*(lam)*(lam))*(5*al - 4))*(D14) + (-14*((lam)*(lam)*(lam))*(5*al - 2))*(D32) + (36*((lam)*(lam)))*(D13) + (36*((lam)*(lam)))*(D31) + (-((al)*(al))*((lam)*(lam)*(lam)*(lam)))*(D60) + (-7*al*((lam)*(lam)*(lam)))*(D50) + (2*al*((lam)*(lam)*(lam)*(lam))*(3*al - 1))*(D51);
        // lambda_alpha_alpha_alpha_alpha_alpha
        o[89] = (-((lam)*(lam)*(lam)*(lam)*(lam))*(al - 1))*(D06) + (((lam)*(lam)*(lam)*(lam)*(lam))*(6*al - 1))*(D51) + (-40*((lam)*(lam)*(lam)*(lam)))*(D32) + (-20*((lam)*(lam)*(lam)*(lam)))*(D14) + (10*((lam)*(lam)*(lam)*(lam)*(lam))*(2*al - 1))*(D33) + (((lam)*(lam)*(lam)*(lam)*(lam))*(6*al - 5))*(D15) + (-4*((lam)*(lam)*(lam)*(lam)))*(D50) + (4*((lam)*(lam)*(lam)*(lam)))*(D05) + (-5*((lam)*(lam)*(lam)*(lam)*(lam))*(3*al - 1))*(D42) + (-5*((lam)*(lam)*(lam)*(lam)*(lam))*(3*al - 2))*(D24) + (20*((lam)*(lam)*(lam)*(lam)))*(D41) + (40*((lam)*(lam)*(lam)*(lam)))*(D23) + (-al*((lam)*(lam)*(lam)*(lam)*(lam)))*(D60);
}


// [[Rcpp::export]]
List enet_dexpected3_cpp(NumericVector y, NumericVector mu,
                         NumericVector lambda, NumericVector alpha,
                         int threads = 1) {
    const std::size_t n = std::max({y.size(), mu.size(), lambda.size(),
                                    alpha.size()});
    std::vector<double> ul, ua, tabD, tabU;
    std::vector<int> id;
    d34e_enet_pairs(lambda, alpha, n, 5, ul, ua, tabD, tabU, id);
    const std::size_t nu = ul.size();
    std::vector<double> comp(nu * 60);
    d7::par_for(nu, threads, d7::kMinCheap, [&](std::size_t u) {
        d34e_enet_comp3(ul[u], ua[u], tabU[u], &tabD[u * 49],
                                &comp[u * 60]);
    });
    return d34e::run(static_cast<int>(n), threads, d7::kMinCostly, 60,
                     [&](std::size_t i, d34e::Row o) {
        const double* c = &comp[id[i] * 60];
        for (int k = 0; k < 60; ++k) o[k] = c[k];
    });
}

// [[Rcpp::export]]
List enet_dexpected4_cpp(NumericVector y, NumericVector mu,
                         NumericVector lambda, NumericVector alpha,
                         int threads = 1) {
    const std::size_t n = std::max({y.size(), mu.size(), lambda.size(),
                                    alpha.size()});
    std::vector<double> ul, ua, tabD, tabU;
    std::vector<int> id;
    d34e_enet_pairs(lambda, alpha, n, 6, ul, ua, tabD, tabU, id);
    const std::size_t nu = ul.size();
    std::vector<double> comp(nu * 90);
    d7::par_for(nu, threads, d7::kMinCheap, [&](std::size_t u) {
        d34e_enet_comp4(ul[u], ua[u], tabU[u], &tabD[u * 49],
                                &comp[u * 90]);
    });
    return d34e::run(static_cast<int>(n), threads, d7::kMinCostly, 90,
                     [&](std::size_t i, d34e::Row o) {
        const double* c = &comp[id[i] * 90];
        for (int k = 0; k < 90; ++k) o[k] = c[k];
    });
}
