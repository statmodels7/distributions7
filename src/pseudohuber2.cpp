#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <map>
#include <string>
#include <vector>
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
// Orders one and two are written out per component. From order three a
// component of D = sqrt(Q) is Faa di Bruno over the set partitions of its
// indices, with (d/dQ)^m sqrt(Q) = c_m Q^(1/2 - m); the partitions of each
// component are enumerated once per order, those with a block of three or
// more mu-derivatives (a zero) dropped and equal ones merged, so that an
// observation costs only products of partials of Q. Each kernel computes
// R and h to its own order and nothing above it.
//
// R and h are read off L_n(nu) = log K_n(sqrt(nu)), n = 1, 2. In t = sqrt(nu)
// the ratios K_n^(j)(t) / K_n(t) are sums of K_m(t) / K_n(t), formed with the
// exponentially scaled Bessel functions so that the e^-t factors cancel; the
// derivatives of log K_n in t follow from those ratios by the moment-cumulant
// relations, and the derivatives in nu from Faa di Bruno through t(nu).

namespace {

// (d/dQ)^m sqrt(Q) = c_m Q^(1/2 - m)
const double SQRT_C[6] = {1.0, 0.5, -0.25, 0.375, -0.9375, 3.28125};

// nu-only quantities to order K: R^(k) and h^(k), k = 0..K (h[0] unused).
struct NuPart {
  double R[6];
  double h[6];
};

void nu_part(double v, int K, NuPart &o) {
  const double t = std::sqrt(v);
  double km[8];
  for (int m = 0; m <= K + 2; m++) km[m] = R::bessel_k(t, (double) m, 2.0);
  o.R[0] = t * km[2] / km[1];
  if (K == 0) return;
  // derivatives of t(nu) = nu^(1/2)
  double tk[6];
  {
    double c = 1.0;
    for (int k = 1; k <= K; k++) {
      c *= (0.5 - (k - 1));
      tk[k] = c * std::pow(v, 0.5 - k);
    }
  }
  double L[3][6];
  for (int n = 1; n <= 2; n++) {
    // rho_j = K_n^(j)(t) / K_n(t)
    double rho[6];
    for (int j = 1; j <= K; j++) {
      double s = 0.0, binom = 1.0;
      for (int i = 0; i <= j; i++) {
        if (i > 0) binom = binom * (j - i + 1) / i;
        s += binom * km[std::abs(n - j + 2 * i)];
      }
      rho[j] = std::pow(-0.5, j) * s / km[n];
    }
    // g_k = (d/dt)^k log K_n(t)
    double g[6];
    const double r1 = rho[1];
    g[1] = r1;
    if (K >= 2) g[2] = rho[2] - r1 * r1;
    if (K >= 3) g[3] = rho[3] - 3.0 * r1 * rho[2] + 2.0 * r1 * r1 * r1;
    if (K >= 4) g[4] = rho[4] - 4.0 * r1 * rho[3] - 3.0 * rho[2] * rho[2] +
                       12.0 * r1 * r1 * rho[2] - 6.0 * std::pow(r1, 4);
    if (K >= 5) g[5] = rho[5] - 5.0 * r1 * rho[4] - 10.0 * rho[2] * rho[3] +
                       20.0 * r1 * r1 * rho[3] + 30.0 * r1 * rho[2] * rho[2] -
                       60.0 * std::pow(r1, 3) * rho[2] + 24.0 * std::pow(r1, 5);
    L[n][1] = g[1] * tk[1];
    if (K >= 2) L[n][2] = g[2] * tk[1] * tk[1] + g[1] * tk[2];
    if (K >= 3) L[n][3] = g[3] * std::pow(tk[1], 3) + 3.0 * g[2] * tk[1] * tk[2] +
                          g[1] * tk[3];
    if (K >= 4) L[n][4] = g[4] * std::pow(tk[1], 4) +
                          6.0 * g[3] * tk[1] * tk[1] * tk[2] +
                          g[2] * (3.0 * tk[2] * tk[2] + 4.0 * tk[1] * tk[3]) +
                          g[1] * tk[4];
    if (K >= 5) L[n][5] = g[5] * std::pow(tk[1], 5) +
                          10.0 * g[4] * std::pow(tk[1], 3) * tk[2] +
                          g[3] * (15.0 * tk[1] * tk[2] * tk[2] +
                                  10.0 * tk[1] * tk[1] * tk[3]) +
                          g[2] * (10.0 * tk[2] * tk[3] + 5.0 * tk[1] * tk[4]) +
                          g[1] * tk[5];
  }
  // (d/dnu)^k log nu = (-1)^(k-1) (k-1)! / nu^k
  double lg[6];
  {
    double f = 1.0;
    for (int k = 1; k <= K; k++) {
      if (k > 1) f *= -(k - 1);
      lg[k] = f / std::pow(v, k);
    }
  }
  double lam[6];
  for (int k = 1; k <= K; k++) lam[k] = 0.5 * lg[k] + L[2][k] - L[1][k];
  // R^(k) = R Y_k(lambda), Y_k the complete Bell polynomial
  o.R[1] = o.R[0] * lam[1];
  if (K >= 2) o.R[2] = o.R[0] * (lam[1] * lam[1] + lam[2]);
  if (K >= 3) o.R[3] = o.R[0] * (std::pow(lam[1], 3) + 3.0 * lam[1] * lam[2] +
                                 lam[3]);
  if (K >= 4) o.R[4] = o.R[0] * (std::pow(lam[1], 4) +
                                 6.0 * lam[1] * lam[1] * lam[2] +
                                 4.0 * lam[1] * lam[3] + 3.0 * lam[2] * lam[2] +
                                 lam[4]);
  if (K >= 5) o.R[5] = o.R[0] * (std::pow(lam[1], 5) +
                                 10.0 * std::pow(lam[1], 3) * lam[2] +
                                 15.0 * lam[1] * lam[2] * lam[2] +
                                 10.0 * lam[1] * lam[1] * lam[3] +
                                 10.0 * lam[2] * lam[3] + 5.0 * lam[1] * lam[4] +
                                 lam[5]);
  for (int k = 1; k <= K; k++)
    o.h[k] = -0.25 * lg[k] + 0.5 * L[2][k] - 1.5 * L[1][k];
}

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
  NuPart np;
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
    const double r = y[i] - m_i;
    const double B[3] = {r * r, -2.0 * r, 2.0};
    double Cj[6];
    const double is = 1.0 / s;
    double sp = is * is;
    for (int j = 0; j <= m; j++) { Cj[j] = cj_num[j] * sp; sp *= is; }
    for (int a = 0; a < 3; a++)
      for (int j = 0; j <= m; j++)
        for (int k = 0; k <= m; k++)
          qt[a * 36 + j * 6 + k] = np.R[k] * B[a] * Cj[j];
    qt[1] += 1.0;   // the nu in Q
    const double Q = qt[0] + v;
    const double D = std::sqrt(Q);
    double dq[6];   // D / Q^b
    dq[0] = D;
    for (int b = 1; b <= m; b++) dq[b] = dq[b - 1] / Q;
    for (int c = 0; c < nc; c++) {
      const Component &cp = C[c];
      double dD = 0.0;
      for (size_t t = 0; t < cp.terms.size(); t++) {
        const Term &tm = cp.terms[t];
        double prod = tm.coef * dq[tm.nb];
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
  NuPart np;
  double last_nu = NA_REAL;
  for (int i = 0; i < n; i++) {
    const double m_i = mu_s ? mu[0] : mu[i];
    const double s = sg_s ? sigma[0] : sigma[i];
    const double v = nu_s ? nu[0] : nu[i];
    if (!(v == last_nu)) { nu_part(v, 1, np); last_nu = v; }
    const double r = y[i] - m_i;
    const double w = r * r / (s * s);
    const double D = std::sqrt(v + np.R[0] * w);
    // l = -D - log sigma + h:  D_x = Q_x / (2 D)
    g_mu[i] = np.R[0] * r / (s * s * D);
    g_sigma[i] = (np.R[0] * w / D - 1.0) / s;
    g_nu[i] = -(1.0 + np.R[1] * w) / (2.0 * D) + np.h[1];
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
  NuPart np;
  double last_nu = NA_REAL;
  for (int i = 0; i < n; i++) {
    const double m_i = mu_s ? mu[0] : mu[i];
    const double s = sg_s ? sigma[0] : sigma[i];
    const double v = nu_s ? nu[0] : nu[i];
    if (!(v == last_nu)) { nu_part(v, 2, np); last_nu = v; }
    const double r = y[i] - m_i;
    const double s2 = s * s, s3 = s2 * s;
    const double R0 = np.R[0], R1 = np.R[1], R2 = np.R[2];
    const double w = r * r / s2;
    const double D = std::sqrt(v + R0 * w);
    const double i2D = 0.5 / D, i4D3 = 0.25 / (D * D * D);
    // partials of Q
    const double Qm = -2.0 * R0 * r / s2, Qs = -2.0 * R0 * w / s,
                 Qn = 1.0 + R1 * w;
    const double Qmm = 2.0 * R0 / s2, Qms = 4.0 * R0 * r / s3,
                 Qmn = -2.0 * R1 * r / s2, Qss = 6.0 * R0 * w / s2,
                 Qsn = -2.0 * R1 * w / s, Qnn = R2 * w;
    // l_xy = -(Q_xy / (2 D) - Q_x Q_y / (4 D^3)) + pure terms
    h_mm[i] = -(Qmm * i2D - Qm * Qm * i4D3);
    h_ss[i] = -(Qss * i2D - Qs * Qs * i4D3) + 1.0 / s2;
    h_nn[i] = -(Qnn * i2D - Qn * Qn * i4D3) + np.h[2];
    h_ms[i] = -(Qms * i2D - Qm * Qs * i4D3);
    h_mn[i] = -(Qmn * i2D - Qm * Qn * i4D3);
    h_sn[i] = -(Qsn * i2D - Qs * Qn * i4D3);
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
    const double t = std::sqrt(v);
    const double k1 = R::bessel_k(t, 1.0, 2.0), k2 = R::bessel_k(t, 2.0, 2.0);
    R[i] = t * k2 / k1;
    // log K_n(t) = log(scaled K_n(t)) - t
    h[i] = -0.25 * std::log(v) + 0.5 * std::log(k2) - 1.5 * std::log(k1) + t;
  }
  return List::create(Named("R") = R, Named("h") = h);
}
