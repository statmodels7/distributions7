#ifndef D7_PT_GAMMA2_H
#define D7_PT_GAMMA2_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <vector>
#include "psi_diff.h"

// Gamma in the mean m and the variance v, with shape a = m^2/v and rate
// m/v. The log-density is written as
//   l = G(a) + a e(m) - log y,   G(a) = a log a - lgamma(a) - a,
//   e(m) = log z - (z - 1),       z = y/m,
// so that every derivative is a sum of products of three kinds of factor:
// G^(r)(a), which are the cancellation-free rests of psi_diff.h
// (G' = log a - psi(a), G'' = 1/a - psi'(a), ...); the partial derivatives
// of the polynomial a(m, v); and e^(k)(m), with e itself formed by psi_Ew2()
// and the higher ones linear in w = (y - m)/m. G(a(m, v)) is differentiated
// by Faa di Bruno over the set partitions of the differentiation slots, and
// a e(m) by the product rule. The direct forms in digamma(a) and log(m/v y)
// lose a digit for every factor of ten in the shape. The expected values
// replace e by its expectation -G'(a) and w by 0: every factor of y enters
// linearly. The vector kernels (gamma.cpp) and the scalar registry
// (d7_ccallable.cpp) call gamma2_component().

namespace d7 {

// One term of Faa di Bruno for G(a(m, v)): G^(nb)(a) times the product over
// the nb blocks of d^(p_b, q_b) a, counted mult times.
struct Gamma2Term {
    int nb;
    int p[4];
    int q[4];
    double mult;
};

// The terms of the derivative of order i in m and j in v, i + j <= 4, from
// the set partitions of the i + j slots (restricted growth strings; slots
// below i are m). Blocks with more than two m slots vanish, a being
// quadratic in m, and partitions with the same blocks are merged.
inline std::vector<Gamma2Term> gamma2_build_terms(int i, int j) {
    const int k = i + j;
    std::map<std::vector<int>, Gamma2Term> acc;
    int rg[4] = {0, 0, 0, 0};
    while (true) {
        int nb = 0;
        for (int s = 0; s < k; ++s) if (rg[s] + 1 > nb) nb = rg[s] + 1;
        std::vector<int> key;
        Gamma2Term t{nb, {0, 0, 0, 0}, {0, 0, 0, 0}, 1.0};
        bool zero = false;
        for (int b = 0; b < nb; ++b) {
            int p = 0, q = 0;
            for (int s = 0; s < k; ++s) {
                if (rg[s] == b) { if (s < i) ++p; else ++q; }
            }
            if (p > 2) zero = true;
            key.push_back(p * 8 + q);
        }
        if (!zero) {
            std::sort(key.begin(), key.end());
            for (int b = 0; b < nb; ++b) { t.p[b] = key[b] / 8; t.q[b] = key[b] % 8; }
            auto it = acc.find(key);
            if (it == acc.end()) acc[key] = t; else it->second.mult += 1.0;
        }
        int t2 = k - 1;
        for (; t2 >= 1; --t2) {
            int mx = 0;
            for (int s = 0; s < t2; ++s) if (rg[s] > mx) mx = rg[s];
            if (rg[t2] <= mx) break;
        }
        if (t2 < 1) break;
        ++rg[t2];
        for (int s = t2 + 1; s < k; ++s) rg[s] = 0;
    }
    std::vector<Gamma2Term> out;
    for (auto& kv : acc) out.push_back(kv.second);
    return out;
}

inline const std::vector<Gamma2Term>& gamma2_terms(int i, int j) {
    static const std::vector<std::vector<std::vector<Gamma2Term>>> tab = [] {
        std::vector<std::vector<std::vector<Gamma2Term>>> t(
            5, std::vector<std::vector<Gamma2Term>>(5));
        for (int i = 0; i <= 4; ++i)
            for (int j = 0; i + j <= 4; ++j)
                if (i + j >= 1) t[i][j] = gamma2_build_terms(i, j);
        return t;
    }();
    return tab[i][j];
}

// Per-observation factors: da[p][q] = d^(p, q) a, de[k] = e^(k)(m), and the
// shape rests Gd[r] = G^(r)(a), r = 1..4.
struct Gamma2Tables {
    double da[3][5];
    double de[5];
    double Gd[5];
};

// The shape rests G^(r)(a), r = 1..order, in Gd[1..order]. For an expected
// value Gd[1] is set to 0: G'(a) and E[e] = -G'(a) multiply the same
// derivative of a in every component, so both are dropped (see
// gamma2_fill_tables()).
inline void gamma2_shape_rests(double a, int order, double* Gd,
                               bool expected = false) {
    Gd[0] = 0.0;
    Gd[1] = expected ? 0.0 : d7::psi_log_rest(a);
    Gd[2] = (order >= 2) ? d7::psi1_rest(a) : 0.0;
    Gd[3] = (order >= 3) ? d7::psi2_rest(a) : 0.0;
    Gd[4] = (order >= 4) ? d7::psi3_rest(a) : 0.0;
}

// da and de up to the given order; Gd must already be filled. With
// expected = true, e is replaced by its expectation -G'(a), which cancels
// G'(a) exactly and is therefore set to 0 together with Gd[1], and w by 0.
inline void gamma2_fill_tables(double y, double m, double v, int order,
                               bool expected, Gamma2Tables& T) {
    const double iv = 1.0 / v;
    const double c[3] = {m * m, 2.0 * m, 2.0};
    for (int p = 0; p < 3; ++p) {
        double f = c[p] * iv;
        for (int q = 0; q <= order; ++q) {
            T.da[p][q] = f;
            f *= -(q + 1.0) * iv;
        }
    }
    double w = 0.0;
    if (expected) {
        T.de[0] = 0.0;
    } else {
        w = (y - m) / m;
        T.de[0] = d7::psi_Ew2(y / m, w);
    }
    // e^(k) = (-1)^k (k - 1)! m^-k ((1 - k) - k w)
    double fk = -1.0 / m;
    for (int k = 1; k <= order; ++k) {
        T.de[k] = fk * ((1.0 - k) - k * w);
        fk *= -static_cast<double>(k) / m;
    }
}

// The derivative of l of order i in m and j in v, from the tables.
inline double gamma2_component(int i, int j, const Gamma2Tables& T) {
    double g = 0.0;
    for (const Gamma2Term& t : gamma2_terms(i, j)) {
        double pr = t.mult * T.Gd[t.nb];
        for (int b = 0; b < t.nb; ++b) pr *= T.da[t.p[b]][t.q[b]];
        g += pr;
    }
    // a e(m) by the product rule in m; e does not depend on v
    double h = 0.0, binom = 1.0;
    for (int kk = 0; kk <= i; ++kk) {
        if (i - kk <= 2) h += binom * T.da[i - kk][j] * T.de[kk];
        binom = binom * (i - kk) / (kk + 1);
    }
    return g + h;
}

inline double gamma2_dexpected_mu_mu_mu(double m, double v, double r1,
                                        double r2) {
    double a = m * m / v;
    double g1 = -r1 - a * r2;
    double iv = 1.0 / v, iv2 = iv * iv;
    return -8.0 * m * g1 * iv2;
}

inline double gamma2_dexpected_sigma2_sigma2_sigma2(double m, double v,
                                                    double r1, double r2) {
    double a = m * m / v;
    double g = -a * r1, g1 = -r1 - a * r2;
    double iv = 1.0 / v, iv2 = iv * iv, iv3 = iv2 * iv;
    return a * (3.0 * g + a * g1) * iv3;
}

// The derivatives of the expected information entries that the mixed
// components of gamma2_dexpected1_cpp hold: dE[mu mu]/dv, dE[v v]/dmu,
// dE[mu v]/dmu and dE[mu v]/dv, as gamma2_dexpected1_cpp computes them
// (r1 = psi1_rest(a), r2 = psi2_rest(a)).
inline double gamma2_dexpected_mu_mu_sigma2(double m, double v, double r1,
                                            double r2) {
    double a = m * m / v;
    double g = -a * r1, g1 = -r1 - a * r2;
    double iv = 1.0 / v, iv2 = iv * iv;
    return (1.0 + 4.0 * g + 4.0 * a * g1) * iv2;
}

inline double gamma2_dexpected_sigma2_sigma2_mu(double m, double v, double r1,
                                                double r2) {
    double a = m * m / v;
    double g = -a * r1, g1 = -r1 - a * r2;
    double iv = 1.0 / v, iv2 = iv * iv, iv3 = iv2 * iv;
    return -2.0 * m * (g + a * g1) * iv3;
}

inline double gamma2_dexpected_mu_sigma2_mu(double m, double v, double r1,
                                            double r2) {
    double a = m * m / v;
    double g = -a * r1, g1 = -r1 - a * r2;
    double iv = 1.0 / v, iv2 = iv * iv;
    return (2.0 * g + 4.0 * a * g1) * iv2;
}

inline double gamma2_dexpected_mu_sigma2_sigma2(double m, double v, double r1,
                                                double r2) {
    double a = m * m / v;
    double g = -a * r1, g1 = -r1 - a * r2;
    double iv = 1.0 / v, iv2 = iv * iv, iv3 = iv2 * iv;
    return -2.0 * m * (2.0 * g + a * g1) * iv3;
}

inline void gamma2_score_curv(int k, double y, const double* th,
                              double* out) {
    double m = th[0], s2 = th[1];
    Gamma2Tables T;
    gamma2_shape_rests(m * m / s2, 2, T.Gd);
    gamma2_fill_tables(y, m, s2, 2, false, T);
    if (k == 0) {
        out[0] = gamma2_component(1, 0, T);
        out[1] = gamma2_component(2, 0, T);
    } else {
        out[0] = gamma2_component(0, 1, T);
        out[1] = gamma2_component(0, 2, T);
    }
}

inline void gamma2_info_dinfo(int k, double y, const double* th,
                              double* out) {
    double m = th[0], s2 = th[1];
    double a = m * m / s2;
    Gamma2Tables T;
    gamma2_shape_rests(a, 2, T.Gd, true);
    gamma2_fill_tables(y, m, s2, 2, true, T);
    double r1 = T.Gd[2], r2 = d7::psi2_rest(a);
    if (k == 0) {
        out[0] = gamma2_component(2, 0, T);
        out[1] = gamma2_dexpected_mu_mu_mu(m, s2, r1, r2);
    } else {
        out[0] = gamma2_component(0, 2, T);
        out[1] = gamma2_dexpected_sigma2_sigma2_sigma2(m, s2, r1, r2);
    }
}

// The routers of the joint registry, one per order: every component at one
// observation, on the parameter scale (mu, sigma2). The components are the
// ones gamma.cpp reads from gamma2_component() on the tables of one
// observation. The Hessian and the expected information fill a full P x P
// matrix and the derivative of the expected information a P x P x P array,
// column-major, out[a + P b + P^2 c] being the derivative of E[l_ab] in
// parameter c.
inline void gamma2_joint_grad(double y, const double* th, double* out) {
    double m = th[0], v = th[1];
    Gamma2Tables T;
    gamma2_shape_rests(m * m / v, 2, T.Gd);
    gamma2_fill_tables(y, m, v, 2, false, T);
    out[0] = gamma2_component(1, 0, T);
    out[1] = gamma2_component(0, 1, T);
}

inline void gamma2_joint_hess(double y, const double* th, double* out) {
    double m = th[0], v = th[1];
    Gamma2Tables T;
    gamma2_shape_rests(m * m / v, 2, T.Gd);
    gamma2_fill_tables(y, m, v, 2, false, T);
    out[0] = gamma2_component(2, 0, T);
    out[3] = gamma2_component(0, 2, T);
    out[1] = out[2] = gamma2_component(1, 1, T);
}

inline void gamma2_joint_info(double y, const double* th, double* out) {
    double m = th[0], v = th[1];
    Gamma2Tables T;
    gamma2_shape_rests(m * m / v, 2, T.Gd, true);
    gamma2_fill_tables(y, m, v, 2, true, T);
    out[0] = gamma2_component(2, 0, T);
    out[3] = gamma2_component(0, 2, T);
    out[1] = out[2] = gamma2_component(1, 1, T);
}

inline void gamma2_joint_dinfo(double y, const double* th, double* out) {
    double m = th[0], v = th[1];
    double a = m * m / v;
    double r1 = d7::psi1_rest(a), r2 = d7::psi2_rest(a);
    for (int j = 0; j < 8; ++j) out[j] = 0.0;
    out[0] = gamma2_dexpected_mu_mu_mu(m, v, r1, r2);
    out[4] = gamma2_dexpected_mu_mu_sigma2(m, v, r1, r2);
    out[1] = out[2] = gamma2_dexpected_mu_sigma2_mu(m, v, r1, r2);
    out[3] = gamma2_dexpected_sigma2_sigma2_mu(m, v, r1, r2);
    out[5] = out[6] = gamma2_dexpected_mu_sigma2_sigma2(m, v, r1, r2);
    out[7] = gamma2_dexpected_sigma2_sigma2_sigma2(m, v, r1, r2);
}

} // namespace d7

#endif
