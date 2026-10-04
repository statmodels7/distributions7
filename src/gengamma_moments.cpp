#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// The central moments of a generalized gamma, on the scale of its mean.
// With k = d/p and h = 1/p, E[Y^j] = E[Y]^j exp(D_j),
//   D_j = lgamma(k + j h) - lgamma(k) - j (lgamma(k + h) - lgamma(k)),
// and X = Y/E[Y] has mu_r = E[(X - 1)^r] = sum_i C(r, i) (-1)^(r-i) e^(D_i).
// Towards the lognormal (k large against h) the terms of that sum agree to
// several orders: written directly the variance lost 2.4e-4 at k = 1e6.
//
// THE SERIES. The cumulant generating function of log X is
//   K(s) = sum_{n>=2} c_n (s^n - s)/n!,  c_n = psi^(n-1)(k) h^n,
// convergent for |s| h < k. With e^K(s) = sum_n a_n s^n, the r-th forward
// difference at zero kills every power below r, so
//   mu_r = r! sum_{n>=r} S(n, r) a_n
// (S the Stirling numbers of the second kind) carries no cancellation of
// the leading orders. The excess kurtosis needs kappa_4 = mu_4 - 3 mu_2^2,
// which cancels again (to 1/k, and to 1/k^2 at p = 2 and p = 4, where the
// coefficient 2(4h - 1)(2h - 1) of 1/k vanishes), and the skewness's mu_3
// cancels to the same extent at p = 3; so the recursion for a_n and the
// sums are carried in double-double, on inputs c_n rounded to double (the
// moments are well conditioned in the c_n; only the arithmetic cancels).
// Measured with stabilita/gg_moments_measure.R against 300 digits.
//
// The series is used for 8h <= k and 4 h^2 <= k, with 60 terms ((4h/k)^60 <
// 1e-18 at the threshold); elsewhere the direct form, with lgamma differences in Stirling's
// form for k >= 10.

namespace {

struct dd {
  double hi, lo;
  dd(double h = 0.0, double l = 0.0) : hi(h), lo(l) {}
};
inline dd two_sum(double a, double b) {
  const double s = a + b, bb = s - a;
  return dd(s, (a - (s - bb)) + (b - bb));
}
inline dd quick(double a, double b) {
  const double s = a + b;
  return dd(s, b - (s - a));
}
inline dd operator+(const dd& a, const dd& b) {
  dd s = two_sum(a.hi, b.hi);
  const dd t = two_sum(a.lo, b.lo);
  s.lo += t.hi;
  s = quick(s.hi, s.lo);
  s.lo += t.lo;
  return quick(s.hi, s.lo);
}
inline dd operator-(const dd& a) { return dd(-a.hi, -a.lo); }
inline dd operator-(const dd& a, const dd& b) { return a + (-b); }
// the exact product by Dekker's split (R builds without -mfma)
inline dd two_prod(double a, double b) {
  const double p = a * b;
  const double ta = 134217729.0 * a, ah = ta - (ta - a), al = a - ah;
  const double tb = 134217729.0 * b, bh = tb - (tb - b), bl = b - bh;
  return dd(p, ((ah * bh - p) + ah * bl + al * bh) + al * bl);
}
inline dd operator*(const dd& a, const dd& b) {
  dd p = two_prod(a.hi, b.hi);
  p.lo += a.hi * b.lo + a.lo * b.hi;
  return quick(p.hi, p.lo);
}
inline dd operator*(double a, const dd& b) { return dd(a) * b; }
inline dd operator/(const dd& a, double b) {
  const double q1 = a.hi / b;
  dd r = a - two_prod(q1, b);
  const double q2 = r.hi / b;
  r = r - two_prod(q2, b);
  const double q3 = r.hi / b;
  return quick(q1, q2) + dd(q3);
}

const int NS = 60;

// the series where it converges fast (4h/k <= 1/2) and the dispersion is
// moderate (h^2/k, about psi'(k) h^2, at most 1/4; measured, at 1 the
// kurtosis lost 2.5e-7 in the series and 8e-13 in the direct form)
inline bool gg_series(double k, double h) { return 8.0 * h <= k && 4.0 * h * h <= k; }

// a_0..a_NS, the coefficients of e^K(s)
void gg_exp_coef(double k, double h, dd* a) {
  double b[NS + 1];
  double fact = 1.0, rk = h / k, rkn = rk, hn = h;
  b[0] = 0.0;
  for (int n = 2; n <= NS; n++) {
    fact *= n;
    rkn *= rk;
    hn *= h;
    // psi^(n-1)(k) = psi^(n-1)(k+1) + (-1)^n (n-1)!/k^n, the two of one
    // sign: the second term is formed as (h/k)^n, which does not overflow
    const double sg = (n % 2) ? -1.0 : 1.0;
    const double cn = R::psigamma(k + 1.0, n - 1.0) * hn + sg * (fact / n) * rkn;
    b[n] = cn / fact;
  }
  dd b1(0.0);
  for (int n = NS; n >= 2; n--) b1 = b1 - dd(b[n]);
  a[0] = dd(1.0);
  for (int m = 1; m <= NS; m++) {
    dd s = (double) 1 * b1 * a[m - 1];
    for (int j = 2; j <= m; j++) s = s + (double) j * (dd(b[j]) * a[m - j]);
    a[m] = s / (double) m;
  }
}

// r! sum_{n>=r} S(n, r) a_n
dd gg_stirling_sum(const dd* a, int r) {
  dd S[NS + 1][5];
  for (int n = 0; n <= NS; n++) for (int j = 0; j <= 4; j++) S[n][j] = dd(0.0);
  S[0][0] = dd(1.0);
  for (int n = 1; n <= NS; n++)
    for (int j = 1; j <= r && j <= n; j++)
      S[n][j] = (double) j * S[n - 1][j] + S[n - 1][j - 1];
  dd s(0.0);
  for (int n = NS; n >= r; n--) s = s + S[n][r] * a[n];
  double rf = 1.0;
  for (int j = 2; j <= r; j++) rf *= j;
  return rf * s;
}

// lgamma(k + x) - lgamma(k), x >= 0
double gg_lgd(double k, double x) {
  if (k < 10.0) return R::lgammafn(k + x) - R::lgammafn(k);
  auto St = [](double z) {
    const double iz = 1.0 / z, iz2 = iz * iz;
    return iz * (1.0 / 12 + iz2 * (-1.0 / 360 + iz2 * (1.0 / 1260 + iz2 * (-1.0 / 1680 +
           iz2 * (1.0 / 1188 + iz2 * (-691.0 / 360360 + iz2 * (1.0 / 156 +
           iz2 * (-3617.0 / 122400))))))));
  };
  return (k - 0.5) * std::log1p(x / k) + x * std::log(k + x) - x + St(k + x) - St(k);
}

// e^(D_j) - 1
double gg_phi(double k, double h, int j) {
  return std::expm1(gg_lgd(k, j * h) - j * gg_lgd(k, h));
}

}  // namespace

// log E[Y/a] = lgamma(k + h) - lgamma(k)
// [[Rcpp::export]]
NumericVector gengamma_logmean_cpp(NumericVector k, NumericVector h) {
  const int n = std::max(k.size(), h.size());
  NumericVector out(n);
  for (int i = 0; i < n; i++) out[i] = gg_lgd(k[i % k.size()], h[i % h.size()]);
  return out;
}

// mu_2 = Var(Y)/E[Y]^2
// [[Rcpp::export]]
NumericVector gengamma_mu2_cpp(NumericVector k, NumericVector h) {
  const int n = std::max(k.size(), h.size());
  NumericVector out(n);
  dd a[NS + 1];
  for (int i = 0; i < n; i++) {
    const double kv = k[i % k.size()], hv = h[i % h.size()];
    if (gg_series(kv, hv)) {
      gg_exp_coef(kv, hv, a);
      out[i] = gg_stirling_sum(a, 2).hi;
    } else {
      out[i] = gg_phi(kv, hv, 2);
    }
  }
  return out;
}

// mu_3 = E[(Y - E[Y])^3]/E[Y]^3
// [[Rcpp::export]]
NumericVector gengamma_mu3_cpp(NumericVector k, NumericVector h) {
  const int n = std::max(k.size(), h.size());
  NumericVector out(n);
  dd a[NS + 1];
  for (int i = 0; i < n; i++) {
    const double kv = k[i % k.size()], hv = h[i % h.size()];
    if (gg_series(kv, hv)) {
      gg_exp_coef(kv, hv, a);
      out[i] = gg_stirling_sum(a, 3).hi;
    } else {
      out[i] = gg_phi(kv, hv, 3) - 3.0 * gg_phi(kv, hv, 2);
    }
  }
  return out;
}

// kappa_4 = mu_4 - 3 mu_2^2, the fourth cumulant of Y/E[Y]
// [[Rcpp::export]]
NumericVector gengamma_kappa4_cpp(NumericVector k, NumericVector h) {
  const int n = std::max(k.size(), h.size());
  NumericVector out(n);
  dd a[NS + 1];
  for (int i = 0; i < n; i++) {
    const double kv = k[i % k.size()], hv = h[i % h.size()];
    if (gg_series(kv, hv)) {
      gg_exp_coef(kv, hv, a);
      const dd m2 = gg_stirling_sum(a, 2);
      out[i] = (gg_stirling_sum(a, 4) - 3.0 * (m2 * m2)).hi;
    } else {
      const double p2 = gg_phi(kv, hv, 2), p3 = gg_phi(kv, hv, 3), p4 = gg_phi(kv, hv, 4);
      out[i] = (p4 - 4.0 * p3 + 6.0 * p2) - 3.0 * p2 * p2;
    }
  }
  return out;
}
