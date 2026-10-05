#include <Rcpp.h>
#include <cmath>
#include <vector>
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
// The quantities of (lambda, alpha) are formed once when both are scalars.

namespace {

struct EnetAt {
    const NumericVector &lam, &al;
    bool both, l_s, a_s;
    d7::EnetPar P0;
    EnetAt(const NumericVector& l, const NumericVector& a)
        : lam(l), al(a), l_s(l.size() == 1), a_s(a.size() == 1) {
        both = l_s && a_s;
        if (both) P0 = d7::enet_par(lam[0], al[0]);
    }
    d7::EnetPar operator()(std::size_t i) const {
        if (both) return P0;
        return d7::enet_par(l_s ? lam[0] : lam[i], a_s ? al[0] : al[i]);
    }
};

}  // namespace

// [[Rcpp::export]]
List enet_gradient_cpp(NumericVector y, NumericVector mu, NumericVector lambda,
                       NumericVector alpha, int threads = 1) {
    int n = y.size();
    NumericVector g_m(n), g_l(n), g_a(n);
    bool m_s = (mu.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin();
    double *o1 = g_m.begin(), *o2 = g_l.begin(), *o3 = g_a.begin();
    EnetAt at(lambda, alpha);
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
    EnetAt at(lambda, alpha);
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
    EnetAt at(lambda, alpha);
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
    EnetAt at(lambda, alpha);
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
