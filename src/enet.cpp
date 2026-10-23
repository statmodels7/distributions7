#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <functional>
#include <string>
#include <unordered_map>
#include "d7_par.h"
#include "pt_enet.h"
using namespace Rcpp;

// The elastic net's score, Hessian, expected information and the first
// derivatives of the expected information, one kernel each, from the
// component functions of pt_enet.h; the off-diagonal components are written
// here. With w = (alpha, 1 - alpha), v = (1, -1), Z'' and Z''' the second
// and third derivatives of log Z in (a, c), d/dlambda = Z-derivative along
// w and d/dalpha = lambda times the Z-derivative along v:
//   l_ml = alpha sgn(z) + (1 - alpha) z,   l_ma = lambda (sgn(z) - z),
//   l_la = (-|z| + z^2/2) - (lambda w'Z''v + Z_a - Z_c),
//   E_ml = E_ma = 0,  E_la = -lambda w'Z''v,
//   E_mm = -(a^2 - 2ac Z_a - 2c^2 Z_c).
// The quantities of (lambda, alpha) are formed once per distinct pair of
// values (enet_unique_pairs()), the pairs being found before the loop.

namespace {

struct PairHash {
    std::size_t operator()(const std::pair<std::uint64_t, std::uint64_t>& k) const {
        return std::hash<std::uint64_t>()(k.first * 0x9E3779B97F4A7C15ULL ^ k.second);
    }
};

// the distinct (lambda, alpha) pairs among the first n elements, compared by
// their bits, and the index of each element's pair
void enet_unique_pairs(const NumericVector& lam, const NumericVector& al,
                       std::size_t n, std::vector<double>& ul,
                       std::vector<double>& ua, std::vector<int>& id) {
    const bool l_s = lam.size() == 1, a_s = al.size() == 1;
    id.assign(n, 0);
    if (l_s && a_s) {
        ul.assign(1, lam[0]);
        ua.assign(1, al[0]);
        return;
    }
    std::unordered_map<std::pair<std::uint64_t, std::uint64_t>, int, PairHash> seen;
    for (std::size_t i = 0; i < n; ++i) {
        const double l = l_s ? lam[0] : lam[i], a = a_s ? al[0] : al[i];
        std::pair<std::uint64_t, std::uint64_t> k;
        std::memcpy(&k.first, &l, sizeof(double));
        std::memcpy(&k.second, &a, sizeof(double));
        auto it = seen.find(k);
        if (it == seen.end()) {
            it = seen.emplace(k, static_cast<int>(ul.size())).first;
            ul.push_back(l);
            ua.push_back(a);
        }
        id[i] = it->second;
    }
}

struct EnetAt {
    std::vector<d7::EnetPar> P;
    std::vector<int> id;
    EnetAt(const NumericVector& l, const NumericVector& a, std::size_t n) {
        std::vector<double> ul, ua;
        enet_unique_pairs(l, a, n, ul, ua, id);
        P.resize(ul.size());
        for (std::size_t u = 0; u < ul.size(); ++u) P[u] = d7::enet_par(ul[u], ua[u]);
    }
    const d7::EnetPar& operator()(std::size_t i) const { return P[id[i]]; }
};

// The derivatives of l = -a|z| - c z^2/2 - log Z(a, c) in (mu, a, c) along
// the index tuple ix (0 = mu, 1 = a, 2 = c) of length m, with s = sgn(z) and
// d[i][j] the partial derivatives of log Z. The data term is quadratic in z
// and linear in each rate, so beyond the second order only d/dmu^2 dc = -1
// carries mu.
double enet_l_ac(int m, const int* ix, double z, double s, double a, double c,
                 const double d[5][5]) {
    int nm = 0, na = 0, nc = 0;
    for (int r = 0; r < m; ++r) {
        if (ix[r] == 0) ++nm;
        else if (ix[r] == 1) ++na;
        else ++nc;
    }
    if (m == 1) {
        if (nm) return a * s + c * z;
        if (na) return -std::fabs(z) - d[1][0];
        return -z * z / 2.0 - d[0][1];
    }
    if (nm == 0) return -d[na][nc];
    if (m == 2) {
        if (nm == 2) return -c;
        return na ? s : z;
    }
    if (m == 3 && nm == 2 && nc == 1) return -1.0;
    return 0.0;
}

// the m-linear form of those derivatives on the directions dir[0..m-1], each
// a vector in (mu, a, c)
double enet_l_form(int m, const double (*dir)[3], double z, double s, double a,
                   double c, const double d[5][5]) {
    int ix[4] = {0, 0, 0, 0};
    double out = 0.0;
    int total = 1;
    for (int r = 0; r < m; ++r) total *= 3;
    for (int t = 0; t < total; ++t) {
        int u = t;
        double w = 1.0;
        for (int r = 0; r < m; ++r) {
            ix[r] = u % 3;
            u /= 3;
            w *= dir[r][ix[r]];
        }
        if (w != 0.0) out += w * enet_l_ac(m, ix, z, s, a, c, d);
    }
    return out;
}

// One derivative of order k in (mu, lambda, alpha), with nm, nl and na the
// counts of each. a = lambda alpha and c = lambda (1 - alpha) are bilinear:
// d(a, c)/dlambda = w = (alpha, 1 - alpha), d(a, c)/dalpha = lambda v,
// d^2(a, c)/dlambda dalpha = v = (1, -1), and nothing above. The partition
// sum keeps the blocks of size one and the pairs {lambda, alpha}: with j such
// pairs, C(nl, j) C(na, j) j! partitions contribute the form of order k - j
// on nm times e_mu, nl - j times w, na - j times lambda v and j times v.
double enet_deriv_lam_al(int nm, int nl, int na, double z, double s,
                         double lam, double al, const double d[5][5]) {
    const double a = lam * al, c = lam * (1.0 - al);
    const double emu[3] = {1.0, 0.0, 0.0}, w[3] = {0.0, al, 1.0 - al},
                 lv[3] = {0.0, lam, -lam}, v[3] = {0.0, 1.0, -1.0};
    static const double fact[5] = {1.0, 1.0, 2.0, 6.0, 24.0};
    double out = 0.0;
    for (int j = 0; j <= std::min(nl, na); ++j) {
        const double coef = fact[nl] / (fact[j] * fact[nl - j]) *
                            fact[na] / (fact[j] * fact[na - j]) * fact[j];
        double dir[4][3];
        int m = 0;
        auto add = [&](const double* e, int times) {
            for (int t = 0; t < times; ++t, ++m)
                for (int q = 0; q < 3; ++q) dir[m][q] = e[q];
        };
        add(emu, nm);
        add(w, nl - j);
        add(lv, na - j);
        add(v, j);
        out += coef * enet_l_form(m, dir, z, s, a, c, d);
    }
    return out;
}

// The derivatives of order k (3 or 4) in (mu, lambda, alpha), keyed as
// deriv_names(). Each is the same for every observation at one parameter
// value except mu_lambda_alpha at order three, sgn(z) - z, so a table is
// formed once per distinct (lambda, alpha) and that one component is added
// per observation.
List enet_deriv_k(NumericVector y, NumericVector mu, NumericVector lambda,
                  NumericVector alpha, int k) {
    const std::size_t n = std::max({y.size(), mu.size(), lambda.size(),
                                    alpha.size()});
    static const char* nmv[3] = {"mu", "lambda", "alpha"};
    // the multi-indices in the order of deriv_names(): nondecreasing tuples
    std::vector<std::array<int, 3>> cnt;
    std::vector<std::string> keys;
    int ix[4];
    std::function<void(int, int)> rec = [&](int pos, int from) {
        if (pos == k) {
            std::array<int, 3> c3 = {0, 0, 0};
            std::string s;
            for (int r = 0; r < k; ++r) {
                ++c3[ix[r]];
                if (r) s += "_";
                s += nmv[ix[r]];
            }
            cnt.push_back(c3);
            keys.push_back(s);
            return;
        }
        for (int v = from; v < 3; ++v) {
            ix[pos] = v;
            rec(pos + 1, v);
        }
    };
    rec(0, 0);
    const int nk = static_cast<int>(cnt.size());

    std::vector<double> ul, ua;
    std::vector<int> id;
    enet_unique_pairs(lambda, alpha, n, ul, ua, id);
    std::vector<double> tab(ul.size() * nk);
    for (std::size_t u = 0; u < ul.size(); ++u) {
        double d[5][5];
        d7::enet_logz_ac(ul[u] * ua[u], ul[u] * (1.0 - ua[u]), k, d);
        for (int c = 0; c < nk; ++c)
            tab[u * nk + c] = enet_deriv_lam_al(cnt[c][0], cnt[c][1], cnt[c][2],
                                                0.0, 0.0, ul[u], ua[u], d);
    }
    // the component that carries the data at order three
    int mla = -1;
    if (k == 3)
        for (int c = 0; c < nk; ++c)
            if (cnt[c][0] == 1 && cnt[c][1] == 1 && cnt[c][2] == 1) mla = c;

    std::vector<NumericVector> cols(nk);
    for (int c = 0; c < nk; ++c) cols[c] = NumericVector(n);
    const bool y_s = y.size() == 1, m_s = mu.size() == 1;
    for (std::size_t i = 0; i < n; ++i) {
        const double* t = &tab[id[i] * nk];
        for (int c = 0; c < nk; ++c) cols[c][i] = t[c];
        if (mla >= 0) {
            const double z = (y_s ? y[0] : y[i]) - (m_s ? mu[0] : mu[i]);
            cols[mla][i] += d7::enet_sign(z) - z;
        }
    }
    List out(nk);
    CharacterVector nm(nk);
    for (int c = 0; c < nk; ++c) {
        out[c] = cols[c];
        nm[c] = keys[c];
    }
    out.attr("names") = nm;
    return out;
}

}  // namespace

// [[Rcpp::export]]
List enet_gradient_cpp(NumericVector y, NumericVector mu, NumericVector lambda,
                       NumericVector alpha, int threads = 1) {
    int n = y.size();
    NumericVector g_m(n), g_l(n), g_a(n);
    bool m_s = (mu.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin();
    double *o1 = g_m.begin(), *o2 = g_l.begin(), *o3 = g_a.begin();
    EnetAt at(lambda, alpha, n);
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        const d7::EnetPar P = at(i);
        double z = yp[i] - (m_s ? mp[0] : mp[i]);
        o1[i] = d7::enet_score_mu(P, z);
        o2[i] = d7::enet_score_lambda(P, z);
        o3[i] = d7::enet_score_alpha(P, z);
    });
    return List::create(Named("mu") = g_m, Named("lambda") = g_l,
                        Named("alpha") = g_a);
}

// [[Rcpp::export]]
List enet_hessian_cpp(NumericVector y, NumericVector mu, NumericVector lambda,
                      NumericVector alpha, int threads = 1) {
    int n = y.size();
    NumericVector h_mm(n), h_ll(n), h_aa(n), h_ml(n), h_ma(n), h_la(n);
    bool m_s = (mu.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin();
    double *o1 = h_mm.begin(), *o2 = h_ll.begin(), *o3 = h_aa.begin(),
           *o4 = h_ml.begin(), *o5 = h_ma.begin(), *o6 = h_la.begin();
    EnetAt at(lambda, alpha, n);
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        const d7::EnetPar P = at(i);
        double z = yp[i] - (m_s ? mp[0] : mp[i]);
        double s = d7::enet_sign(z), al = P.al, lam = P.lam;
        o1[i] = d7::enet_hess_mu_mu(P);
        o2[i] = d7::enet_hess_lambda_lambda(P);
        o3[i] = d7::enet_hess_alpha_alpha(P);
        o4[i] = al * s + (1.0 - al) * z;
        o5[i] = lam * (s - z);
        o6[i] = (-std::fabs(z) + z * z / 2.0) -
            (lam * (P.zaa * al + P.zac * (1.0 - 2.0 * al) - P.zcc * (1.0 - al)) +
             P.za - P.zc);
    });
    return List::create(Named("mu_mu") = h_mm, Named("lambda_lambda") = h_ll,
                        Named("alpha_alpha") = h_aa, Named("mu_lambda") = h_ml,
                        Named("mu_alpha") = h_ma, Named("lambda_alpha") = h_la);
}

// [[Rcpp::export]]
List enet_expected_hessian_cpp(NumericVector y, NumericVector mu,
                               NumericVector lambda, NumericVector alpha,
                               int threads = 1) {
    int n = y.size();
    NumericVector h_mm(n), h_ll(n), h_aa(n), h_ml(n), h_ma(n), h_la(n);
    double *o1 = h_mm.begin(), *o2 = h_ll.begin(), *o3 = h_aa.begin(),
           *o6 = h_la.begin();
    EnetAt at(lambda, alpha, n);
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        const d7::EnetPar P = at(i);
        double al = P.al;
        o1[i] = d7::enet_expected_mu_mu(P);
        o2[i] = d7::enet_expected_lambda_lambda(P);
        o3[i] = d7::enet_expected_alpha_alpha(P);
        o6[i] = -P.lam * (P.zaa * al + P.zac * (1.0 - 2.0 * al) -
                          P.zcc * (1.0 - al));
    });
    return List::create(Named("mu_mu") = h_mm, Named("lambda_lambda") = h_ll,
                        Named("alpha_alpha") = h_aa, Named("mu_lambda") = h_ml,
                        Named("mu_alpha") = h_ma, Named("lambda_alpha") = h_la);
}

// The first derivatives of the expected information, keyed as
// dexpected_names(): every derivative in mu is zero, and so is every
// derivative of E_ml and E_ma.
// [[Rcpp::export]]
List enet_dexpected1_cpp(NumericVector y, NumericVector mu,
                         NumericVector lambda, NumericVector alpha,
                         int threads = 1) {
    int n = y.size();
    const char* keys[18] = {
        "mu_mu_mu", "mu_mu_lambda", "mu_mu_alpha",
        "lambda_lambda_mu", "lambda_lambda_lambda", "lambda_lambda_alpha",
        "alpha_alpha_mu", "alpha_alpha_lambda", "alpha_alpha_alpha",
        "mu_lambda_mu", "mu_lambda_lambda", "mu_lambda_alpha",
        "mu_alpha_mu", "mu_alpha_lambda", "mu_alpha_alpha",
        "lambda_alpha_mu", "lambda_alpha_lambda", "lambda_alpha_alpha"};
    std::vector<NumericVector> v(18);
    std::vector<double*> o(18);
    for (int k = 0; k < 18; ++k) {
        v[k] = NumericVector(n);
        o[k] = v[k].begin();
    }
    EnetAt at(lambda, alpha, n);
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        const d7::EnetPar P = at(i);
        const d7::EnetT3 T = d7::enet_t3(P);
        double lam = P.lam, al = P.al, bl = 1.0 - al, a = P.a, c = P.c;
        // Z'' along (w, v) and (v, v)
        double Hwv = P.zaa * al + P.zac * (bl - al) - P.zcc * bl;
        double Hvv = P.zaa - 2.0 * P.zac + P.zcc;
        // Z''' along (w, w, v), (w, v, w) = (w, w, v), (v, v, w), (w, v, v)
        double Twwv = T.aaa * al * al + T.aac * (2.0 * al * bl - al * al) +
                      T.acc * (bl * bl - 2.0 * al * bl) - T.ccc * bl * bl;
        double Tvvw = T.aaa * al + T.aac * (bl - 2.0 * al) +
                      T.acc * (al - 2.0 * bl) + T.ccc * bl;
        // E_mm = -(a^2 - 2ac Z_a - 2c^2 Z_c), differentiated in a and in c
        double dEa = -(2.0 * a - 2.0 * c * P.za - 2.0 * a * c * P.zaa -
                       2.0 * c * c * P.zac);
        double dEc = -(-2.0 * a * P.za - 2.0 * a * c * P.zac - 4.0 * c * P.zc -
                       2.0 * c * c * P.zcc);
        o[1][i] = al * dEa + bl * dEc;
        o[2][i] = lam * (dEa - dEc);
        o[4][i] = d7::enet_dexpected_lambda_lambda_lambda(P, T);
        o[5][i] = -(2.0 * Hwv + lam * Twwv);
        o[7][i] = -2.0 * lam * Hvv - lam * lam * Tvvw;
        o[8][i] = d7::enet_dexpected_alpha_alpha_alpha(P, T);
        o[16][i] = -Hwv - lam * Twwv;
        o[17][i] = -lam * (Hvv + lam * Tvvw);
    });
    List out(18);
    CharacterVector nm(18);
    for (int k = 0; k < 18; ++k) {
        out[k] = v[k];
        nm[k] = keys[k];
    }
    out.attr("names") = nm;
    return out;
}

// The second derivatives of the expected information, keyed as
// d2expected_names(): the pair of the information first and the pair of the
// derivatives second, each in the order mu_mu, lambda_lambda, alpha_alpha,
// mu_lambda, mu_alpha, lambda_alpha. Only E_mm, E_ll, E_aa and E_la depend on
// the parameters, and only on (lambda, alpha), so the nine derivatives of each
// (enet_d2expected_block()) are formed once per distinct pair and the other
// 24 keys are zero.
// [[Rcpp::export]]
List enet_dexpected2_cpp(NumericVector y, NumericVector mu,
                         NumericVector lambda, NumericVector alpha,
                         int threads = 1) {
    const int n = y.size();
    static const char* hn[6] = {"mu_mu", "lambda_lambda", "alpha_alpha",
                                "mu_lambda", "mu_alpha", "lambda_alpha"};
    // the row of the block for each pair of the information (-1: zero), and
    // the column for each pair of the derivatives (-1: zero)
    static const int rowOf[6] = {0, 1, 2, -1, -1, 3};
    static const int colOf[6] = {-1, 0, 1, -1, -1, 2};
    std::vector<double> ul, ua;
    std::vector<int> id;
    enet_unique_pairs(lambda, alpha, n, ul, ua, id);
    const std::size_t nu = ul.size();
    std::vector<double> tab(nu * 12);
    for (std::size_t u = 0; u < nu; ++u) {
        double blk[4][3];
        d7::enet_d2expected_block(ul[u], ua[u], blk);
        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 3; ++c) tab[u * 12 + r * 3 + c] = blk[r][c];
    }
    std::vector<NumericVector> v(36);
    std::vector<double*> o(36);
    for (int k = 0; k < 36; ++k) {
        v[k] = NumericVector(n);
        o[k] = v[k].begin();
    }
    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        const double* t = &tab[id[i] * 12];
        for (int ab = 0; ab < 6; ++ab) {
            if (rowOf[ab] < 0) continue;
            for (int cd = 0; cd < 6; ++cd)
                if (colOf[cd] >= 0)
                    o[ab * 6 + cd][i] = t[rowOf[ab] * 3 + colOf[cd]];
        }
    });
    List out(36);
    CharacterVector nm(36);
    for (int ab = 0; ab < 6; ++ab)
        for (int cd = 0; cd < 6; ++cd) {
            out[ab * 6 + cd] = v[ab * 6 + cd];
            nm[ab * 6 + cd] = std::string(hn[ab]) + "_" + hn[cd];
        }
    out.attr("names") = nm;
    return out;
}

// L^(k)(x), k = 1..order, of the logarithm of the Mills ratio
// [[Rcpp::export]]
List enet_logM_derivs_cpp(NumericVector x, int order) {
    const int n = x.size();
    List out(order);
    std::vector<NumericVector> cols;
    for (int k = 0; k < order; k++) cols.push_back(NumericVector(n));
    for (int i = 0; i < n; i++) {
        double L[6];
        d7::enet_logM_derivs(x[i], order, L);
        for (int k = 0; k < order; k++) cols[k][i] = L[k + 1];
    }
    for (int k = 0; k < order; k++) out[k] = cols[k];
    return out;
}

// the partial derivatives of log Z in (a, c) to the given order, as a list
// of orders, each keyed as deriv_names(c("a", "c"), k)
// [[Rcpp::export]]
List enet_logz_ac_cpp(NumericVector a, NumericVector c, int order) {
    const int n = std::max(a.size(), c.size());
    const int na = a.size(), nc = c.size();
    List out(order);
    std::vector<std::vector<NumericVector>> cols(order + 1);
    for (int k = 1; k <= order; k++)
        for (int j = 0; j <= k; j++) cols[k].push_back(NumericVector(n));
    for (int i = 0; i < n; i++) {
        double d[5][5];
        d7::enet_logz_ac(a[i % na], c[i % nc], order, d);
        for (int k = 1; k <= order; k++)
            for (int j = 0; j <= k; j++) cols[k][j][i] = d[k - j][j];
    }
    for (int k = 1; k <= order; k++) {
        List ok(k + 1);
        CharacterVector nm(k + 1);
        for (int j = 0; j <= k; j++) {
            std::string s;
            for (int m = 0; m < k; m++) {
                if (m) s += "_";
                s += (m < k - j) ? "a" : "c";
            }
            ok[j] = cols[k][j];
            nm[j] = s;
        }
        ok.attr("names") = nm;
        out[k - 1] = ok;
    }
    return out;
}

// the third derivatives of the log-density, keyed as
// deriv_names(c("mu", "lambda", "alpha"), 3)
// [[Rcpp::export]]
List enet_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector lambda,
                     NumericVector alpha) {
    return enet_deriv_k(y, mu, lambda, alpha, 3);
}

// the fourth derivatives of the log-density, keyed as
// deriv_names(c("mu", "lambda", "alpha"), 4)
// [[Rcpp::export]]
List enet_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector lambda,
                     NumericVector alpha) {
    return enet_deriv_k(y, mu, lambda, alpha, 4);
}
