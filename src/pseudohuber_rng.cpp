#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "pt_sqrt.h"
using namespace Rcpp;

// Random generation for the pseudo-Huber (symmetric hyperbolic) family as a
// normal variance mixture: with the scale sigma1 of pseudohuber_distrib(),
//   Y = mu + sqrt(W) Z,   W = sigma1^2 sqrt(nu) X,   X ~ GIG(lambda = 1, omega),
// omega = sqrt(nu), where GIG(lambda, omega) has the density proportional to
//   x^(lambda - 1) exp(-omega (x + 1/x) / 2),   x > 0.
// E[W] = sigma1^2 sqrt(nu) K_2(omega) / K_1(omega), the family's variance.
//
// X is drawn by the ratio-of-uniforms method with the mode shifted to the
// origin (Hormann and Leydold, 2014, the region lambda >= 1). For lambda = 1
// the mode is 1, and with h(x) = exp(-omega (x - 1)^2 / (2 x)) = g(x) / g(1)
// the acceptance region is { (u, v): 0 < v <= sqrt(h(u / v + 1)) }, enclosed
// in [u-, u+] x [0, 1] with u = (x - 1) sqrt(h(x)) at the two extrema x-, x+
// of that function. They are the roots in (0, 1) and (1, inf) of
//   x^3 - (1 + 4/omega) x^2 - x + 1 = 0,
// the condition d/dx log((x - 1) sqrt(h(x))) = 0 multiplied out. For omega
// below 1e-8 the GIG is the exponential limit, W = 2 sigma1^2 E.

namespace {

struct RouBox {
  double um, up;
};

double cubic(double x, double a) { return ((x + a) * x - 1.0) * x + 1.0; }
double cubic_d(double x, double a) { return (3.0 * x + 2.0 * a) * x - 1.0; }

RouBox rou_box(double omega) {
  const double a = -(1.0 + 4.0 / omega), b = -1.0, c = 1.0;
  // depressed cubic t^3 + p t + q, x = t - a/3
  const double p = b - a * a / 3.0;
  const double q = 2.0 * a * a * a / 27.0 - a * b / 3.0 + c;
  double arg = -q / 2.0 * d7::sqrt_cr(-27.0 / (p * p * p));
  arg = std::max(-1.0, std::min(1.0, arg));
  const double phi = std::acos(arg);
  const double r = d7::sqrt_cr(-4.0 * p / 3.0);
  double xs[3];
  for (int k = 0; k < 3; k++)
    xs[k] = r * std::cos((phi + 2.0 * M_PI * k) / 3.0) - a / 3.0;
  double xm = NA_REAL, xp = NA_REAL;
  for (int k = 0; k < 3; k++) {
    if (xs[k] > 0.0 && xs[k] < 1.0) xm = xs[k];
    if (xs[k] > 1.0) xp = xs[k];
  }
  // two Newton steps on the cubic sharpen each root
  for (int it = 0; it < 2; it++) {
    xm -= cubic(xm, a) / cubic_d(xm, a);
    xp -= cubic(xp, a) / cubic_d(xp, a);
  }
  RouBox box;
  box.um = (xm - 1.0) * std::exp(-omega * (xm - 1.0) * (xm - 1.0) / (4.0 * xm));
  box.up = (xp - 1.0) * std::exp(-omega * (xp - 1.0) * (xp - 1.0) / (4.0 * xp));
  return box;
}

double gig1_draw(double omega, const RouBox &box) {
  while (true) {
    const double u = box.um + (box.up - box.um) * unif_rand();
    const double v = unif_rand();
    const double x = u / v + 1.0;
    if (x <= 0.0) continue;
    // accept when v^2 <= h(x), on the log scale
    if (2.0 * std::log(v) <= -omega * (x - 1.0) * (x - 1.0) / (2.0 * x)) return x;
  }
}

}  // namespace

// Draws Y for scale sigma1, shape nu and location mu, each recycled to n.
// [[Rcpp::export]]
NumericVector pseudohuber_rng_cpp(int n, NumericVector mu, NumericVector sigma1,
                                  NumericVector nu) {
  RNGScope scope;
  NumericVector out(n);
  const int nm = mu.size(), ns = sigma1.size(), nn = nu.size();
  double last_omega = NA_REAL;
  RouBox box = {0.0, 0.0};
  for (int i = 0; i < n; i++) {
    const double m = mu[i % nm], s1 = sigma1[i % ns], v = nu[i % nn];
    const double omega = d7::sqrt_cr(v);
    double w;
    if (omega < 1e-8) {
      w = 2.0 * s1 * s1 * exp_rand();
    } else {
      if (!(omega == last_omega)) {
        box = rou_box(omega);
        last_omega = omega;
      }
      w = s1 * s1 * omega * gig1_draw(omega, box);
    }
    out[i] = m + d7::sqrt_cr(w) * norm_rand();
  }
  return out;
}

// The bounding box, for the tests.
// [[Rcpp::export]]
NumericVector pseudohuber_rou_box_cpp(double omega) {
  RouBox b = rou_box(omega);
  return NumericVector::create(b.um, b.up);
}
