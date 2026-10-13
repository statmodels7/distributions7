#include <Rcpp.h>
#include <cmath>
#include <string>
#include <vector>
#include "pt_fdb2.h"
using namespace Rcpp;

// The observed derivatives of the lognormal in the mean m and the variance v
// of Y, orders three to five. With S = log(1 + v/m^2) and C = log(y/m) the
// log-density is
//   l = h1(S) + h2(S) C^2 - C/2 - log y - log(2 pi)/2,
//   h1(S) = -log(S)/2 - S/8,   h2(S) = -1/(2S),
// which is the closed form -log(S)/2 - (C + S/2)^2/(2S) with the square
// expanded. Each term is differentiated as a composition: S = log1p(r) with
// r = v/m^2 a monomial, h1 and h2 of S by Faa di Bruno (pt_fdb2.h), and C^2
// by Leibniz from C^(k) = (-1)^k (k-1)!/m^k. Every factor stays of the size
// of the quantity it stands for when v/m^2 is small, where the expanded
// closed forms in (C, S) cancel: at m = 100 and a coefficient of variation
// of 1 per cent they lose eleven digits at order five. Orders one and two
// keep the closed forms of lognormal2.cpp, which lose at most three digits
// there and cost a fraction of this.

namespace {

// Everything that depends on the parameters only, to order K.
struct Ln2Par {
    double dH1[6][6];
    double dH2[6][6];
    double Cm[6];   // C^(k) for k >= 1; Cm[0] unused
    double m;
};

void ln2_par(double m, double v, int K, Ln2Par& P) {
    P.m = m;
    const double im = 1.0 / m, r = v * im * im;
    // r = v m^-2 and its partials: d^p m^-2 = (-1)^p (p+1)! m^-(2+p)
    double dr[6][6] = {{0.0}};
    double fm = im * im;
    for (int p = 0; p <= K; ++p) {
        dr[p][0] = v * fm;
        dr[p][1] = fm;
        fm *= -(p + 2.0) * im;
    }
    // S = log1p(r): s^(k) = (-1)^(k-1) (k-1)! / (1 + r)^k
    double sd[6] = {0.0};
    {
        const double i1 = 1.0 / (1.0 + r);
        double f = i1;
        for (int k = 1; k <= K; ++k) { sd[k] = f; f *= -k * i1; }
    }
    double dS[6][6] = {{0.0}};
    const double S = std::log1p(r);
    dS[0][0] = S;
    for (int p = 0; p <= K; ++p)
        for (int q = 0; p + q <= K; ++q)
            if (p + q >= 1) dS[p][q] = d7::fdb2_compose(p, q, sd, dr);
    // h1(S) = -log(S)/2 - S/8 and h2(S) = -1/(2S)
    double h1d[6] = {0.0}, h2d[6] = {0.0};
    {
        const double iS = 1.0 / S;
        double f = -0.5 * iS;              // -(1/2) (-1)^(k-1) (k-1)! S^-k
        double g = 0.5 * iS * iS;          // -(1/2) (-1)^k k! S^-(k+1)
        for (int k = 1; k <= K; ++k) {
            h1d[k] = f;
            h2d[k] = g;
            f *= -k * iS;
            g *= -(k + 1.0) * iS;
        }
        h1d[1] -= 0.125;
    }
    // h1 enters only at order K; h2 enters at every order up to K, through
    // the product with C^2
    for (int i = 0; i <= K; ++i)
        for (int j = 0; i + j <= K; ++j) {
            if (i + j == 0) { P.dH1[0][0] = 0.0; P.dH2[0][0] = -0.5 / S; continue; }
            P.dH1[i][j] = (i + j == K) ? d7::fdb2_compose(i, j, h1d, dS) : 0.0;
            P.dH2[i][j] = d7::fdb2_compose(i, j, h2d, dS);
        }
    // C^(k) = (-1)^k (k-1)! / m^k
    double c = -im;
    P.Cm[0] = 0.0;
    for (int k = 1; k <= K; ++k) { P.Cm[k] = c; c *= -k * im; }
}

// The (i, j) component at one observation, from C = log(y/m).
inline double ln2_component(int i, int j, const Ln2Par& P, const double* C2d,
                            double C) {
    double s = P.dH1[i][j];
    double binom = 1.0;
    for (int kk = 0; kk <= i; ++kk) {
        s += binom * P.dH2[i - kk][j] * C2d[kk];
        binom = binom * (i - kk) / (kk + 1);
    }
    if (j == 0) s -= 0.5 * (i == 0 ? C : P.Cm[i]);
    return s;
}

// (C^2)^(k) by Leibniz, k = 0..K
inline void ln2_c2(const Ln2Par& P, double C, int K, double* C2d) {
    double Cd[6];
    Cd[0] = C;
    for (int k = 1; k <= K; ++k) Cd[k] = P.Cm[k];
    for (int k = 0; k <= K; ++k) {
        double s = 0.0, binom = 1.0;
        for (int l = 0; l <= k; ++l) {
            s += binom * Cd[l] * Cd[k - l];
            binom = binom * (k - l) / (l + 1);
        }
        C2d[k] = s;
    }
}

// The components (i, j) of one order, named as the R side expects.
List ln2_order(NumericVector y, NumericVector mean, NumericVector var, int K,
               const std::vector<std::pair<int, int>>& idx) {
    const int n = y.size();
    const int n_m = mean.size(), n_v = var.size();
    d7::fdb2_terms(1, 0);
    std::vector<NumericVector> out(idx.size());
    for (auto& o : out) o = NumericVector(n);
    Ln2Par P;
    const bool scalar = n_m == 1 && n_v == 1;
    if (scalar) ln2_par(mean[0], var[0], K, P);
    double C2d[6];
    for (int i = 0; i < n; ++i) {
        if (!scalar) ln2_par(mean[i % n_m], var[i % n_v], K, P);
        // log1p near y = m; elsewhere the difference of logarithms, since
        // (y - m)/m rounds to -1 once y/m is below eps and log1p(-1) = -Inf
        const double C = std::fabs(y[i] - P.m) < 0.5 * P.m
                             ? std::log1p((y[i] - P.m) / P.m)
                             : std::log(y[i]) - std::log(P.m);
        ln2_c2(P, C, K, C2d);
        for (std::size_t c = 0; c < idx.size(); ++c)
            out[c][i] = ln2_component(idx[c].first, idx[c].second, P, C2d, C);
    }
    List res(idx.size());
    CharacterVector nm(idx.size());
    for (std::size_t c = 0; c < idx.size(); ++c) {
        res[c] = out[c];
        std::string s;
        for (int a = 0; a < idx[c].first; ++a) s += (s.empty() ? "" : "_") + std::string("mean");
        for (int b = 0; b < idx[c].second; ++b) s += (s.empty() ? "" : "_") + std::string("var");
        nm[c] = s;
    }
    res.attr("names") = nm;
    return res;
}

// (i, j) in the order of deriv_names(): i from K down to 0
std::vector<std::pair<int, int>> ln2_lex(int K) {
    std::vector<std::pair<int, int>> v;
    for (int i = K; i >= 0; --i) v.push_back({i, K - i});
    return v;
}

} // namespace

// order 3 of log f in (mean, var)
// [[Rcpp::export]]
List lognormal2_deriv3_cpp(NumericVector y, NumericVector mean, NumericVector var) {
    return ln2_order(y, mean, var, 3, ln2_lex(3));
}

// order 4 of log f in (mean, var)
// [[Rcpp::export]]
List lognormal2_deriv4_cpp(NumericVector y, NumericVector mean, NumericVector var) {
    return ln2_order(y, mean, var, 4, ln2_lex(4));
}

// order 5 of log f in (mean, var)
// [[Rcpp::export]]
List lognormal2_deriv5_cpp(NumericVector y, NumericVector mean, NumericVector var) {
    return ln2_order(y, mean, var, 5, ln2_lex(5));
}
