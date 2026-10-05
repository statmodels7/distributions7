#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <map>
#include <string>
#include <vector>
#include "pt_pseudohuber2.h"
#include "pt_sqrt.h"
using namespace Rcpp;

// The pseudo-Huber family in its standard-deviation parametrization
// (mu, sigma, nu), sigma being the standard deviation of Y. With
//   R(nu) = sqrt(nu) K_2(sqrt(nu)) / K_1(sqrt(nu)),   r = y - mu,
//   Q = nu + R(nu) r^2 / sigma^2,   D = sqrt(Q),
// the log-density is
//   l = -D - log 2 - log sigma + h(nu),
//   h(nu) = -log(nu)/4 + log K_2(sqrt(nu))/2 - 3 log K_1(sqrt(nu))/2.
//
// Q is a product of one-variable functions, so every partial of Q is closed:
//   d^i_mu d^j_sigma d^k_nu Q = [i=j=0, k=1] + R^(k)(nu) B_i C_j,
//   B_0 = r^2, B_1 = -2 r, B_2 = 2, B_i = 0 for i >= 3,
//   C_j = (-1)^j (j+1)! sigma^-(j+2).
// Orders one and two are written out per component in pt_pseudohuber2.h,
// with the diagonal of order three the quadrature reads. From order three a
// component of D = sqrt(Q) is Faa di Bruno over the set partitions of its
// indices, with (d/dQ)^m sqrt(Q) = c_m Q^(1/2 - m); the partitions of each
// component are enumerated once per order, those with a block of three or
// more mu-derivatives (a zero) dropped and equal ones merged, so that an
// observation costs only products of partials of Q. Each kernel computes
// R and h to its own order and nothing above it.
//
// R and h (nu_part(), in pt_pseudohuber2.h) are read off
// L_n(nu) = log K_n(sqrt(nu)), n = 1, 2. In t = sqrt(nu)
// the ratios K_n^(j)(t) / K_n(t) are sums of K_m(t) / K_n(t), formed with the
// exponentially scaled Bessel functions so that the e^-t factors cancel; the
// derivatives of log K_n in t follow from those ratios by the moment-cumulant
// relations, and the derivatives in nu from Faa di Bruno through t(nu).

namespace {

using d7::NuPart;
using d7::nu_part;

// (d/dQ)^m sqrt(Q) = c_m Q^(1/2 - m)
const double SQRT_C[6] = {1.0, 0.5, -0.25, 0.375, -0.9375, 3.28125};

// ---- orders three to five: precompiled Faa di Bruno sums ------------------

// A block of a partition is coded i * 36 + j * 6 + k (i mu-, j sigma- and
// k nu-derivatives of Q).
struct Term {
  double coef;               // multiplicity times c_{blocks}
  int nb;                    // number of blocks
  int code[5];
};
struct Component {
  std::vector<Term> terms;
  int nmu, nsig, nnu;
};

void set_partitions(int m, int pos, std::vector<std::vector<int> > &cur,
                    std::vector<std::vector<std::vector<int> > > &out) {
  if (pos == m) { out.push_back(cur); return; }
  for (size_t b = 0; b < cur.size(); b++) {
    cur[b].push_back(pos);
    set_partitions(m, pos + 1, cur, out);
    cur[b].pop_back();
  }
  cur.push_back(std::vector<int>(1, pos));
  set_partitions(m, pos + 1, cur, out);
  cur.pop_back();
}

// the components of order m in lexicographic order of the parameter
// positions, each with its merged, zero-free partition terms
const std::vector<Component> &compiled(int m) {
  static std::vector<Component> cache[6];
  if (!cache[m].empty()) return cache[m];
  std::vector<std::vector<std::vector<int> > > parts;
  std::vector<std::vector<int> > cur;
  set_partitions(m, 0, cur, parts);
  std::vector<int> idx(m, 0);
  while (true) {
    Component comp;
    comp.nmu = comp.nsig = comp.nnu = 0;
    for (int e = 0; e < m; e++) {
      if (idx[e] == 0) comp.nmu++; else if (idx[e] == 1) comp.nsig++; else comp.nnu++;
    }
    std::map<std::vector<int>, double> merged;
    for (size_t p = 0; p < parts.size(); p++) {
      std::vector<int> codes;
      bool zero = false;
      for (size_t b = 0; b < parts[p].size(); b++) {
        int c[3] = {0, 0, 0};
        for (size_t e = 0; e < parts[p][b].size(); e++) c[idx[parts[p][b][e]]]++;
        if (c[0] >= 3) { zero = true; break; }
        codes.push_back(c[0] * 36 + c[1] * 6 + c[2]);
      }
      if (zero) continue;
      std::sort(codes.begin(), codes.end());
      merged[codes] += SQRT_C[codes.size()];
    }
    for (std::map<std::vector<int>, double>::iterator it = merged.begin();
         it != merged.end(); ++it) {
      Term t;
      t.coef = it->second;
      t.nb = (int) it->first.size();
      for (int b = 0; b < t.nb; b++) t.code[b] = it->first[b];
      comp.terms.push_back(t);
    }
    cache[m].push_back(comp);
    int k = m - 1;
    while (k >= 0 && idx[k] == 2) k--;
    if (k < 0) break;
    idx[k]++;
    for (int j = k + 1; j < m; j++) idx[j] = idx[k];
  }
  return cache[m];
}

std::string comp_name(const std::vector<int> &idx) {
  const char *nm[3] = {"mu", "sigma", "nu"};
  std::string key;
  for (size_t e = 0; e < idx.size(); e++) {
    if (e) key += "_";
    key += nm[idx[e]];
  }
  return key;
}

List high_order(int m, NumericVector y, NumericVector mu, NumericVector sigma,
                NumericVector nu) {
  const int n = y.size();
  const bool mu_s = mu.size() == 1, sg_s = sigma.size() == 1,
             nu_s = nu.size() == 1;
  const std::vector<Component> &C = compiled(m);
  const int nc = C.size();
  std::vector<NumericVector> res;
  for (int c = 0; c < nc; c++) res.push_back(NumericVector(n));
  NuPart np = {{0}, {0}};
  double last_nu = NA_REAL;
  // (j+1)! (-1)^j for C_j
  double cj_num[6];
  {
    double f = 1.0;
    for (int j = 0; j <= m; j++) {
      f *= (j + 1);
      cj_num[j] = ((j % 2) ? -1.0 : 1.0) * f;
    }
  }
  // (d/dsigma)^m (-log sigma) = (-1)^m (m-1)! / sigma^m
  double fm = 1.0;
  for (int k = 2; k < m; k++) fm *= k;
  fm *= (m % 2) ? -1.0 : 1.0;
  double qt[3 * 36];
  for (int i = 0; i < n; i++) {
    const double m_i = mu_s ? mu[0] : mu[i];
    const double s = sg_s ? sigma[0] : sigma[i];
    const double v = nu_s ? nu[0] : nu[i];
    if (!(v == last_nu)) { nu_part(v, m, np); last_nu = v; }
    // every partial of Q enters divided by Q, so that a product of up to
    // five partials, each of the order of z^2, never forms: with
    // z = (y - mu)/sigma, B_a C_j / Q is z^2/Q sigma^-j, -2 z/Q sigma^-(j+1)
    // and 2/Q sigma^-(j+2), and z^2/Q = 1/(R + nu/z^2) is bounded. D is
    // ph2_z()'s, and the product is D times the ratios.
    const double z = (y[i] - m_i) / s;
    const d7::Ph2Z P = d7::ph2_z(z, v, np.R[0]);
    double zq2, iQ;   // z^2/Q and 1/Q
    if (std::fabs(z) > 1.0) {
      zq2 = 1.0 / (np.R[0] + v / (z * z));
      iQ = zq2 / (z * z);
    } else {
      iQ = 1.0 / (v + np.R[0] * (z * z));
      zq2 = z * z * iQ;
    }
    const double Bq[3] = {zq2, -2.0 * z * iQ, 2.0 * iQ};
    double Cj[6];
    const double is = 1.0 / s;
    for (int j = 0; j <= m; j++) Cj[j] = cj_num[j] * std::pow(is, j);
    for (int a = 0; a < 3; a++)
      for (int j = 0; j <= m; j++)
        for (int k = 0; k <= m; k++)
          qt[a * 36 + j * 6 + k] = np.R[k] * Bq[a] * Cj[j] * std::pow(is, a);
    qt[1] += iQ;   // the nu in Q
    for (int c = 0; c < nc; c++) {
      const Component &cp = C[c];
      double dD = 0.0;
      for (size_t t = 0; t < cp.terms.size(); t++) {
        const Term &tm = cp.terms[t];
        double prod = tm.coef * P.D;
        for (int b = 0; b < tm.nb; b++) prod *= qt[tm.code[b]];
        dD += prod;
      }
      double out = -dD;
      if (cp.nsig == m) out += fm * std::pow(is, m);
      if (cp.nnu == m) out += np.h[m];
      res[c][i] = out;
    }
  }
  // names from the same enumeration
  List out(nc);
  CharacterVector names(nc);
  std::vector<int> idx(m, 0);
  for (int c = 0; c < nc; c++) {
    out[c] = res[c];
    names[c] = comp_name(idx);
    int k = m - 1;
    while (k >= 0 && idx[k] == 2) k--;
    if (k < 0) break;
    idx[k]++;
    for (int j = k + 1; j < m; j++) idx[j] = idx[k];
  }
  out.attr("names") = names;
  return out;
}

}  // namespace

// [[Rcpp::export]]
List pseudohuber2_gradient_cpp(NumericVector y, NumericVector mu,
                               NumericVector sigma, NumericVector nu) {
  const int n = y.size();
  const bool mu_s = mu.size() == 1, sg_s = sigma.size() == 1,
             nu_s = nu.size() == 1;
  NumericVector g_mu(n), g_sigma(n), g_nu(n);
  NuPart np = {{0}, {0}};
  double last_nu = NA_REAL;
  for (int i = 0; i < n; i++) {
    const double m_i = mu_s ? mu[0] : mu[i];
    const double s = sg_s ? sigma[0] : sigma[i];
    const double v = nu_s ? nu[0] : nu[i];
    if (!(v == last_nu)) { nu_part(v, 1, np); last_nu = v; }
    const double z = (y[i] - m_i) / s;
    const d7::Ph2Z P = d7::ph2_z(z, v, np.R[0]);
    g_mu[i] = d7::pseudohuber2_score_mu(s, P, np);
    g_sigma[i] = d7::pseudohuber2_score_sigma(z, s, P, np);
    g_nu[i] = d7::pseudohuber2_score_nu(z, P, np);
  }
  return List::create(Named("mu") = g_mu, Named("sigma") = g_sigma,
                      Named("nu") = g_nu);
}

// [[Rcpp::export]]
List pseudohuber2_hessian_cpp(NumericVector y, NumericVector mu,
                              NumericVector sigma, NumericVector nu) {
  const int n = y.size();
  const bool mu_s = mu.size() == 1, sg_s = sigma.size() == 1,
             nu_s = nu.size() == 1;
  NumericVector h_mm(n), h_ss(n), h_nn(n), h_ms(n), h_mn(n), h_sn(n);
  NuPart np = {{0}, {0}};
  double last_nu = NA_REAL;
  for (int i = 0; i < n; i++) {
    const double m_i = mu_s ? mu[0] : mu[i];
    const double s = sg_s ? sigma[0] : sigma[i];
    const double v = nu_s ? nu[0] : nu[i];
    if (!(v == last_nu)) { nu_part(v, 2, np); last_nu = v; }
    const double z = (y[i] - m_i) / s;
    const double R0 = np.R[0], R1 = np.R[1];
    const d7::Ph2Z P = d7::ph2_z(z, v, R0);
    const double zD = P.zD, zz = z * zD, s2 = s * s;
    double qn, qnD;
    d7::ph2_qn(z, P, np, &qn, &qnD);
    // l_xy = -(Q_xy / (2 D) - (Q_x / D)(Q_y / D) / (4 D)) + pure terms,
    // with Q_mu / D = -2 R zD / sigma and Q_sigma / D = -2 R z zD / sigma
    h_mm[i] = d7::pseudohuber2_hess_mu_mu(s, P, np, v);
    h_ss[i] = d7::pseudohuber2_hess_sigma_sigma(z, s, P, np, v);
    h_nn[i] = d7::pseudohuber2_hess_nu_nu(z, P, np);
    const double vD2 = v / (P.D * P.D);
    // 2 - R zD^2 = 1 + nu/D^2 and 1 - R zD^2 / 2 = (1 + nu/D^2)/2 exactly
    h_ms[i] = -R0 * zD * (1.0 + vD2) / s2;
    h_mn[i] = zD * (0.5 * R1 * (1.0 + vD2) - 0.5 * R0 / (P.D * P.D)) / s;
    h_sn[i] = zz * (0.5 * R1 * (1.0 + vD2) - 0.5 * R0 / (P.D * P.D)) / s;
  }
  return List::create(Named("mu_mu") = h_mm, Named("sigma_sigma") = h_ss,
                      Named("nu_nu") = h_nn, Named("mu_sigma") = h_ms,
                      Named("mu_nu") = h_mn, Named("sigma_nu") = h_sn);
}

// [[Rcpp::export]]
List pseudohuber2_deriv3_cpp(NumericVector y, NumericVector mu,
                             NumericVector sigma, NumericVector nu) {
  return high_order(3, y, mu, sigma, nu);
}

// [[Rcpp::export]]
List pseudohuber2_deriv4_cpp(NumericVector y, NumericVector mu,
                             NumericVector sigma, NumericVector nu) {
  return high_order(4, y, mu, sigma, nu);
}

// [[Rcpp::export]]
List pseudohuber2_deriv5_cpp(NumericVector y, NumericVector mu,
                             NumericVector sigma, NumericVector nu) {
  return high_order(5, y, mu, sigma, nu);
}

// R(nu) and the log normalizing term h(nu), for the density in R.
// [[Rcpp::export]]
List pseudohuber2_nu_terms_cpp(NumericVector nu) {
  const int n = nu.size();
  NumericVector R(n), h(n);
  for (int i = 0; i < n; i++) {
    const double v = nu[i];
    const d7::Ph2Terms T = d7::pseudohuber2_terms(v);
    R[i] = T.R;
    h[i] = T.h;
  }
  return List::create(Named("R") = R, Named("h") = h);
}
