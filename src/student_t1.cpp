#include <Rcpp.h>
#include <cmath>
#include "d7_par.h"
#include "pt_student_t1.h"
using namespace Rcpp;

// The Student t in its location mu, scale s and degrees of freedom nu > 0.
// With z = (y - mu)/s, k = nu, q = z^2/k and t = 1/(1 + q), the
// log-density is
//   c(nu) - log s - log(pi)/2 - (k + 1)/2 log(1 + q),
//   c(nu) = lgamma((nu+1)/2) - lgamma(nu/2) - log(nu)/2.
// Every component below is a closed form derived offline with sympy
// (stabilita/gen_student_t2.py, T_FAMILY=t1), one exported function per
// order and per surface, and is checked numerically in the tests.
//
// Every derivative in nu vanishes as nu grows and is a difference of terms
// agreeing to leading order. The data part is cancelled symbolically and
// written in q, z, 1/k and t, which neither cancel nor overflow at the nu
// the link can produce; the score's logarithm enters as
// D(q) = q/(1+q) - log1p(q). Each quantity of nu alone that carries a
// polygamma, t1_V*, has two branches: the direct form below nu = 20, and
// above it the asymptotic series in h = 1/nu, exact to h^40, from
// Stirling's series. The expectations are closed, t being Beta(nu/2, 1/2)
// under the model.

// D(q) and the quantities of nu alone, t1_*, are in pt_student_t1.h.

// order 1 of log f in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t1_gradient_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu(n);
  NumericVector o_sigma(n);
  NumericVector o_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, V0; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double V0 = t1_V0(v);
    return Prm{m, s, ik, V0};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double V0 = P.V0; (void) V0;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    (void) t; (void) q;
    const double DQ = d7::student_t1_DQ(z, ik);
    o_mu[i] = d7::student_t1_score_mu(z, s, ik);
    o_sigma[i] = d7::student_t1_score_sigma(z, s, ik);
    o_nu[i] = d7::student_t1_score_nu(z, ik, V0, DQ);
  });
  return List::create(Named("mu") = o_mu, Named("sigma") = o_sigma, Named("nu") = o_nu);
}

// order 2 of log f in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t1_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu(n);
  NumericVector o_sigma_sigma(n);
  NumericVector o_nu_nu(n);
  NumericVector o_mu_sigma(n);
  NumericVector o_mu_nu(n);
  NumericVector o_sigma_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, V1, w1, w6; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double V1 = t1_V1(v);
    const double w1 = std::pow(s, -2);
    const double w6 = std::pow(ik, 2);
    return Prm{m, s, ik, V1, w1, w6};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double V1 = P.V1; (void) V1;
    const double w1 = P.w1; (void) w1;
    const double w6 = P.w6; (void) w6;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    // u = q/(1 + q) in [0, 1], and zt = z t, where q itself may overflow
    const double u = (q < 1.0) ? q * t : 1.0 / (1.0 + 1.0 / q);
    const double zt = (std::fabs(z) <= 1.0) ? z * t : 1.0 / (1.0 / z + z * ik);
    (void) t; (void) q; (void) u; (void) zt;
    const double w0 = ik*t - u;
    const double w2 = t*w1;
    const double w7 = -w0/s;
    o_mu_mu[i] = d7::student_t1_hess_mu_mu(z, s, ik);
    o_sigma_sigma[i] = d7::student_t1_hess_sigma_sigma(z, s, ik);
    o_nu_nu[i] = d7::student_t1_hess_nu_nu(z, ik, V1);
    o_mu_sigma[i] = -2*w2*zt*(ik + 1);
    o_mu_nu[i] = ik*w7*zt;
    o_sigma_nu[i] = u*w7;
  });
  return List::create(Named("mu_mu") = o_mu_mu, Named("sigma_sigma") = o_sigma_sigma, Named("nu_nu") = o_nu_nu, Named("mu_sigma") = o_mu_sigma, Named("mu_nu") = o_mu_nu, Named("sigma_nu") = o_sigma_nu);
}

// order 2, expected
// [[Rcpp::export]]
List student_t1_expected_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu(n);
  NumericVector o_sigma_sigma(n);
  NumericVector o_nu_nu(n);
  NumericVector o_mu_sigma(n);
  NumericVector o_mu_nu(n);
  NumericVector o_sigma_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, V2, w0, w1; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double V2 = t1_V2(v);
    const double w0 = 1/(std::pow(s, 2)*(3*ik + 1));
    const double w1 = std::pow(ik, 2);
    return Prm{m, s, ik, V2, w0, w1};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double V2 = P.V2; (void) V2;
    const double w0 = P.w0; (void) w0;
    const double w1 = P.w1; (void) w1;
    o_mu_mu[i] = d7::student_t1_expected_mu_mu(s, ik);
    o_sigma_sigma[i] = d7::student_t1_expected_sigma_sigma(s, ik);
    o_nu_nu[i] = d7::student_t1_expected_nu_nu(V2);
    o_mu_sigma[i] = 0;
    o_mu_nu[i] = 0;
    o_sigma_nu[i] = 2*w1/(s*(4*ik + 3*w1 + 1));
  });
  return List::create(Named("mu_mu") = o_mu_mu, Named("sigma_sigma") = o_sigma_sigma, Named("nu_nu") = o_nu_nu, Named("mu_sigma") = o_mu_sigma, Named("mu_nu") = o_mu_nu, Named("sigma_nu") = o_sigma_nu);
}

// order 3 of log f in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t1_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu_mu(n);
  NumericVector o_mu_mu_sigma(n);
  NumericVector o_mu_mu_nu(n);
  NumericVector o_mu_sigma_sigma(n);
  NumericVector o_mu_sigma_nu(n);
  NumericVector o_mu_nu_nu(n);
  NumericVector o_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_nu(n);
  NumericVector o_sigma_nu_nu(n);
  NumericVector o_nu_nu_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, V3, w5, w11, w14, w20; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double V3 = t1_V3(v);
    const double w5 = 2/std::pow(s, 3);
    const double w11 = std::pow(s, -2);
    const double w14 = std::pow(ik, 2);
    const double w20 = std::pow(ik, 3);
    return Prm{m, s, ik, V3, w5, w11, w14, w20};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double V3 = P.V3; (void) V3;
    const double w5 = P.w5; (void) w5;
    const double w11 = P.w11; (void) w11;
    const double w14 = P.w14; (void) w14;
    const double w20 = P.w20; (void) w20;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    // u = q/(1 + q) in [0, 1], and zt = z t, where q itself may overflow
    const double u = (q < 1.0) ? q * t : 1.0 / (1.0 + 1.0 / q);
    const double zt = (std::fabs(z) <= 1.0) ? z * t : 1.0 / (1.0 / z + z * ik);
    (void) t; (void) q; (void) u; (void) zt;
    const double w0 = 3*t;
    const double w1 = -u;
    const double w2 = ik*u;
    const double w3 = ik*w0 + w0 + w1 - w2;
    const double w4 = ik*t;
    const double w6 = w5*zt;
    const double w7 = std::pow(t, 2);
    const double w8 = 3*u;
    const double w9 = std::pow(u, 2);
    const double w10 = u*w0;
    const double w12 = w11*w4;
    const double w13 = 2*zt;
    const double w15 = t*w14;
    const double w16 = (w1 + w4)/s;
    const double w17 = w7*w8;
    const double w18 = z*zt;
    const double w19 = t*w2;
    const double w21 = std::pow(u, 3);
    o_mu_mu_mu[i] = -w3*w4*w6;
    o_mu_mu_sigma[i] = w5*w7*(t - 3*w2 + w4 - w8);
    o_mu_mu_nu[i] = w12*(ik*w7 - w0*w2 - w10 + w9);
    o_mu_sigma_sigma[i] = t*w3*w6;
    o_mu_sigma_nu[i] = w12*w13*(ik*t - 2*u - w2);
    o_mu_nu_nu[i] = w13*w15*w16;
    o_sigma_sigma_sigma[i] = w5*(-std::pow(t, 3) + w10*w18 + w17 + 6*w18*w7 + w18*w9);
    o_sigma_sigma_nu[i] = u*w11*(3*ik*w7 - 5*t*u - w19 - w9);
    o_sigma_nu_nu[i] = 2*w16*w19;
    o_nu_nu_nu[i] = V3 + w0*w20*w9 - 1.0/2.0*w14*w21 - 3.0/2.0*w15*w9 + w17*w20 + w20*w21;
  });
  return List::create(Named("mu_mu_mu") = o_mu_mu_mu, Named("mu_mu_sigma") = o_mu_mu_sigma, Named("mu_mu_nu") = o_mu_mu_nu, Named("mu_sigma_sigma") = o_mu_sigma_sigma, Named("mu_sigma_nu") = o_mu_sigma_nu, Named("mu_nu_nu") = o_mu_nu_nu, Named("sigma_sigma_sigma") = o_sigma_sigma_sigma, Named("sigma_sigma_nu") = o_sigma_sigma_nu, Named("sigma_nu_nu") = o_sigma_nu_nu, Named("nu_nu_nu") = o_nu_nu_nu);
}

// order 3, expected
// [[Rcpp::export]]
List student_t1_deriv3_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu_mu(n);
  NumericVector o_mu_mu_sigma(n);
  NumericVector o_mu_mu_nu(n);
  NumericVector o_mu_sigma_sigma(n);
  NumericVector o_mu_sigma_nu(n);
  NumericVector o_mu_nu_nu(n);
  NumericVector o_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_nu(n);
  NumericVector o_sigma_nu_nu(n);
  NumericVector o_nu_nu_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, V4, w0, w1, w2, w3, w4, w5, w6; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double V4 = t1_V4(v);
    const double w0 = std::pow(ik, 2);
    const double w1 = 2*w0;
    const double w2 = 1.0/(8*ik + 15*w0 + 1);
    const double w3 = 2*w2/std::pow(s, 3);
    const double w4 = ik - 1;
    const double w5 = w2/std::pow(s, 2);
    const double w6 = std::pow(ik, 3);
    return Prm{m, s, ik, V4, w0, w1, w2, w3, w4, w5, w6};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double V4 = P.V4; (void) V4;
    const double w0 = P.w0; (void) w0;
    const double w1 = P.w1; (void) w1;
    const double w2 = P.w2; (void) w2;
    const double w3 = P.w3; (void) w3;
    const double w4 = P.w4; (void) w4;
    const double w5 = P.w5; (void) w5;
    const double w6 = P.w6; (void) w6;
    o_mu_mu_mu[i] = 0;
    o_mu_mu_sigma[i] = w3*(3*ik + w1 + 1);
    o_mu_mu_nu[i] = w1*w4*w5;
    o_mu_sigma_sigma[i] = 0;
    o_mu_sigma_nu[i] = 0;
    o_mu_nu_nu[i] = 0;
    o_sigma_sigma_sigma[i] = w3*(13*ik + 5);
    o_sigma_sigma_nu[i] = -12*w0*w5;
    o_sigma_nu_nu[i] = 4*w4*w6/(s*(9*ik + 23*w0 + 15*w6 + 1));
    o_nu_nu_nu[i] = V4;
  });
  return List::create(Named("mu_mu_mu") = o_mu_mu_mu, Named("mu_mu_sigma") = o_mu_mu_sigma, Named("mu_mu_nu") = o_mu_mu_nu, Named("mu_sigma_sigma") = o_mu_sigma_sigma, Named("mu_sigma_nu") = o_mu_sigma_nu, Named("mu_nu_nu") = o_mu_nu_nu, Named("sigma_sigma_sigma") = o_sigma_sigma_sigma, Named("sigma_sigma_nu") = o_sigma_sigma_nu, Named("sigma_nu_nu") = o_sigma_nu_nu, Named("nu_nu_nu") = o_nu_nu_nu);
}

// order 4 of log f in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t1_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu_mu_mu(n);
  NumericVector o_mu_mu_mu_sigma(n);
  NumericVector o_mu_mu_mu_nu(n);
  NumericVector o_mu_mu_sigma_sigma(n);
  NumericVector o_mu_mu_sigma_nu(n);
  NumericVector o_mu_mu_nu_nu(n);
  NumericVector o_mu_sigma_sigma_sigma(n);
  NumericVector o_mu_sigma_sigma_nu(n);
  NumericVector o_mu_sigma_nu_nu(n);
  NumericVector o_mu_nu_nu_nu(n);
  NumericVector o_sigma_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_sigma_nu(n);
  NumericVector o_sigma_sigma_nu_nu(n);
  NumericVector o_sigma_nu_nu_nu(n);
  NumericVector o_nu_nu_nu_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, V5, w0, w16, w17, w19, w22, w23, w26, w36, w37; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double V5 = t1_V5(v);
    const double w0 = std::pow(s, -4);
    const double w16 = std::pow(ik, 2);
    const double w17 = 2/std::pow(s, 3);
    const double w19 = 6*w0;
    const double w22 = std::pow(s, -2);
    const double w23 = w16*w22;
    const double w26 = std::pow(ik, 3);
    const double w36 = std::pow(ik, 4);
    const double w37 = 12*w36;
    return Prm{m, s, ik, V5, w0, w16, w17, w19, w22, w23, w26, w36, w37};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double V5 = P.V5; (void) V5;
    const double w0 = P.w0; (void) w0;
    const double w16 = P.w16; (void) w16;
    const double w17 = P.w17; (void) w17;
    const double w19 = P.w19; (void) w19;
    const double w22 = P.w22; (void) w22;
    const double w23 = P.w23; (void) w23;
    const double w26 = P.w26; (void) w26;
    const double w36 = P.w36; (void) w36;
    const double w37 = P.w37; (void) w37;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    // u = q/(1 + q) in [0, 1], and zt = z t, where q itself may overflow
    const double u = (q < 1.0) ? q * t : 1.0 / (1.0 + 1.0 / q);
    const double zt = (std::fabs(z) <= 1.0) ? z * t : 1.0 / (1.0 / z + z * ik);
    (void) t; (void) q; (void) u; (void) zt;
    const double w1 = ik*t;
    const double w2 = u*w1;
    const double w3 = -6*w2;
    const double w4 = t*u;
    const double w5 = -6*w4;
    const double w6 = std::pow(u, 2);
    const double w7 = ik*w6;
    const double w8 = std::pow(t, 2);
    const double w9 = ik*w8;
    const double w10 = w6 + w9;
    const double w11 = w10 + w3 + w5 + w7 + w8;
    const double w12 = 6*w9;
    const double w13 = -u + w1;
    const double w14 = -ik*u + t + w13;
    const double w15 = 24*w0*zt;
    const double w18 = w17*zt;
    const double w20 = 6*w6;
    const double w21 = -8*ik*t*u;
    const double w24 = 3*w9;
    const double w25 = 4*t;
    const double w27 = w26*w8;
    const double w28 = -6*w13/s;
    const double w29 = std::pow(t, 3);
    const double w30 = 6*w29;
    const double w31 = std::pow(u, 3);
    const double w32 = z*zt;
    const double w33 = w25*w6;
    const double w34 = u*w8;
    const double w35 = std::pow(u, 4);
    o_mu_mu_mu_mu[i] = w0*w11*w12;
    o_mu_mu_mu_sigma[i] = w14*w15*w9;
    o_mu_mu_mu_nu[i] = t*w16*w18*(w12 + w3 - 8*w4 + w6 + 3*w8);
    o_mu_mu_sigma_sigma[i] = -w11*w19*w8;
    o_mu_mu_sigma_nu[i] = w17*w9*(-w20 - w21 - w5 - 3*w7 - w9);
    o_mu_mu_nu_nu[i] = 2*w23*w8*(5*ik*t*u + 3*t*u - 3*w6 - w9);
    o_mu_sigma_sigma_sigma[i] = -w14*w15*w8;
    o_mu_sigma_sigma_nu[i] = w1*w18*(10*t*u - w21 - w24 - 2*w6 - w7);
    o_mu_sigma_nu_nu[i] = w23*w25*zt*(2*ik*t*u + 2*t*u - w10);
    o_mu_nu_nu_nu[i] = w27*w28*zt;
    o_sigma_sigma_sigma_sigma[i] = w19*(std::pow(t, 4) - u*w30 - 10*w29*w32 - w31*w32 - w32*w33 - 5*w32*w34 + w6*w8);
    o_sigma_sigma_sigma_nu[i] = u*w17*(-ik*w30 + u*w12 + w31 + w33 + 15*w34);
    o_sigma_sigma_nu_nu[i] = 2*w2*w22*(3*ik*t*u + 5*t*u - w24 - w6);
    o_sigma_nu_nu_nu[i] = w16*w28*w34;
    o_nu_nu_nu_nu[i] = V5 - t*w31*w37 - u*w29*w37 + w20*w27 + w25*w26*w31 + w26*w35 - 3*w35*w36 - 18*w36*w6*w8;
  });
  return List::create(Named("mu_mu_mu_mu") = o_mu_mu_mu_mu, Named("mu_mu_mu_sigma") = o_mu_mu_mu_sigma, Named("mu_mu_mu_nu") = o_mu_mu_mu_nu, Named("mu_mu_sigma_sigma") = o_mu_mu_sigma_sigma, Named("mu_mu_sigma_nu") = o_mu_mu_sigma_nu, Named("mu_mu_nu_nu") = o_mu_mu_nu_nu, Named("mu_sigma_sigma_sigma") = o_mu_sigma_sigma_sigma, Named("mu_sigma_sigma_nu") = o_mu_sigma_sigma_nu, Named("mu_sigma_nu_nu") = o_mu_sigma_nu_nu, Named("mu_nu_nu_nu") = o_mu_nu_nu_nu, Named("sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma, Named("sigma_sigma_sigma_nu") = o_sigma_sigma_sigma_nu, Named("sigma_sigma_nu_nu") = o_sigma_sigma_nu_nu, Named("sigma_nu_nu_nu") = o_sigma_nu_nu_nu, Named("nu_nu_nu_nu") = o_nu_nu_nu_nu);
}

// order 4, expected
// [[Rcpp::export]]
List student_t1_deriv4_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu_mu_mu(n);
  NumericVector o_mu_mu_mu_sigma(n);
  NumericVector o_mu_mu_mu_nu(n);
  NumericVector o_mu_mu_sigma_sigma(n);
  NumericVector o_mu_mu_sigma_nu(n);
  NumericVector o_mu_mu_nu_nu(n);
  NumericVector o_mu_sigma_sigma_sigma(n);
  NumericVector o_mu_sigma_sigma_nu(n);
  NumericVector o_mu_sigma_nu_nu(n);
  NumericVector o_mu_nu_nu_nu(n);
  NumericVector o_sigma_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_sigma_nu(n);
  NumericVector o_sigma_sigma_nu_nu(n);
  NumericVector o_sigma_nu_nu_nu(n);
  NumericVector o_nu_nu_nu_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, V6, w0, w1, w2, w3, w4, w5, w6, w7, w8, w9; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double V6 = t1_V6(v);
    const double w0 = std::pow(ik, 2);
    const double w1 = 2*w0;
    const double w2 = 1/(std::pow(s, 4)*(12*ik + 35*w0 + 1));
    const double w3 = 6*w2*(3*ik + w1 + 1);
    const double w4 = std::pow(ik, 3);
    const double w5 = 1.0/(15*ik + 71*w0 + 105*w4 + 1);
    const double w6 = w5/std::pow(s, 3);
    const double w7 = 4*w0 - 1;
    const double w8 = w4*w5/std::pow(s, 2);
    const double w9 = std::pow(ik, 4);
    return Prm{m, s, ik, V6, w0, w1, w2, w3, w4, w5, w6, w7, w8, w9};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double V6 = P.V6; (void) V6;
    const double w0 = P.w0; (void) w0;
    const double w1 = P.w1; (void) w1;
    const double w2 = P.w2; (void) w2;
    const double w3 = P.w3; (void) w3;
    const double w4 = P.w4; (void) w4;
    const double w5 = P.w5; (void) w5;
    const double w6 = P.w6; (void) w6;
    const double w7 = P.w7; (void) w7;
    const double w8 = P.w8; (void) w8;
    const double w9 = P.w9; (void) w9;
    o_mu_mu_mu_mu[i] = ik*w3;
    o_mu_mu_mu_sigma[i] = 0;
    o_mu_mu_mu_nu[i] = 0;
    o_mu_mu_sigma_sigma[i] = -w3;
    o_mu_mu_sigma_nu[i] = w1*w6*(9*ik - w1 + 5);
    o_mu_mu_nu_nu[i] = -4*w7*w8;
    o_mu_sigma_sigma_sigma[i] = 0;
    o_mu_sigma_sigma_nu[i] = 0;
    o_mu_sigma_nu_nu[i] = 0;
    o_mu_nu_nu_nu[i] = 0;
    o_sigma_sigma_sigma_sigma[i] = -18*w2*(11*ik + 3);
    o_sigma_sigma_sigma_nu[i] = 6*w0*w6*(31*ik + 13);
    o_sigma_sigma_nu_nu[i] = -12*w8*(ik - 2);
    o_sigma_nu_nu_nu[i] = -12*w7*w9/(s*(16*ik + 86*w0 + 176*w4 + 105*w9 + 1));
    o_nu_nu_nu_nu[i] = V6;
  });
  return List::create(Named("mu_mu_mu_mu") = o_mu_mu_mu_mu, Named("mu_mu_mu_sigma") = o_mu_mu_mu_sigma, Named("mu_mu_mu_nu") = o_mu_mu_mu_nu, Named("mu_mu_sigma_sigma") = o_mu_mu_sigma_sigma, Named("mu_mu_sigma_nu") = o_mu_mu_sigma_nu, Named("mu_mu_nu_nu") = o_mu_mu_nu_nu, Named("mu_sigma_sigma_sigma") = o_mu_sigma_sigma_sigma, Named("mu_sigma_sigma_nu") = o_mu_sigma_sigma_nu, Named("mu_sigma_nu_nu") = o_mu_sigma_nu_nu, Named("mu_nu_nu_nu") = o_mu_nu_nu_nu, Named("sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma, Named("sigma_sigma_sigma_nu") = o_sigma_sigma_sigma_nu, Named("sigma_sigma_nu_nu") = o_sigma_sigma_nu_nu, Named("sigma_nu_nu_nu") = o_sigma_nu_nu_nu, Named("nu_nu_nu_nu") = o_nu_nu_nu_nu);
}

// order 5 of log f in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t1_deriv5_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu_mu_mu_mu(n);
  NumericVector o_mu_mu_mu_mu_sigma(n);
  NumericVector o_mu_mu_mu_mu_nu(n);
  NumericVector o_mu_mu_mu_sigma_sigma(n);
  NumericVector o_mu_mu_mu_sigma_nu(n);
  NumericVector o_mu_mu_mu_nu_nu(n);
  NumericVector o_mu_mu_sigma_sigma_sigma(n);
  NumericVector o_mu_mu_sigma_sigma_nu(n);
  NumericVector o_mu_mu_sigma_nu_nu(n);
  NumericVector o_mu_mu_nu_nu_nu(n);
  NumericVector o_mu_sigma_sigma_sigma_sigma(n);
  NumericVector o_mu_sigma_sigma_sigma_nu(n);
  NumericVector o_mu_sigma_sigma_nu_nu(n);
  NumericVector o_mu_sigma_nu_nu_nu(n);
  NumericVector o_mu_nu_nu_nu_nu(n);
  NumericVector o_sigma_sigma_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_sigma_sigma_nu(n);
  NumericVector o_sigma_sigma_sigma_nu_nu(n);
  NumericVector o_sigma_sigma_nu_nu_nu(n);
  NumericVector o_sigma_nu_nu_nu_nu(n);
  NumericVector o_nu_nu_nu_nu_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, V7, w1, w11, w25, w26, w33, w35, w36, w43, w44, w45, w48, w57, w59, w60, w63; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double V7 = t1_V7(v);
    const double w1 = std::pow(ik, 2);
    const double w11 = 24/std::pow(s, 5);
    const double w25 = std::pow(s, -4);
    const double w26 = 6*w25;
    const double w33 = std::pow(ik, 3);
    const double w35 = std::pow(s, -3);
    const double w36 = 12*w35;
    const double w43 = 4*w35;
    const double w44 = std::pow(s, -2);
    const double w45 = w33*w44;
    const double w48 = std::pow(ik, 4);
    const double w57 = std::pow(ik, 5);
    const double w59 = 60*w57;
    const double w60 = 30*w48;
    const double w63 = 120*w57;
    return Prm{m, s, ik, V7, w1, w11, w25, w26, w33, w35, w36, w43, w44, w45, w48, w57, w59, w60, w63};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double V7 = P.V7; (void) V7;
    const double w1 = P.w1; (void) w1;
    const double w11 = P.w11; (void) w11;
    const double w25 = P.w25; (void) w25;
    const double w26 = P.w26; (void) w26;
    const double w33 = P.w33; (void) w33;
    const double w35 = P.w35; (void) w35;
    const double w36 = P.w36; (void) w36;
    const double w43 = P.w43; (void) w43;
    const double w44 = P.w44; (void) w44;
    const double w45 = P.w45; (void) w45;
    const double w48 = P.w48; (void) w48;
    const double w57 = P.w57; (void) w57;
    const double w59 = P.w59; (void) w59;
    const double w60 = P.w60; (void) w60;
    const double w63 = P.w63; (void) w63;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    // u = q/(1 + q) in [0, 1], and zt = z t, where q itself may overflow
    const double u = (q < 1.0) ? q * t : 1.0 / (1.0 + 1.0 / q);
    const double zt = (std::fabs(z) <= 1.0) ? z * t : 1.0 / (1.0 / z + z * ik);
    (void) t; (void) q; (void) u; (void) zt;
    const double w0 = std::pow(t, 2);
    const double w2 = w0*w1;
    const double w3 = std::pow(u, 2);
    const double w4 = ik*w3;
    const double w5 = 5*w0;
    const double w6 = t*u;
    const double w7 = ik*t;
    const double w8 = u*w7;
    const double w9 = -10*w6 - 10*w8;
    const double w10 = ik*w5 + w3 + w4 + w5 + w9;
    const double w12 = w11*zt;
    const double w13 = w10*w12;
    const double w14 = 5*w3;
    const double w15 = ik*w0;
    const double w16 = w14 + w15;
    const double w17 = w0 + w16 + 5*w4 + w9;
    const double w18 = std::pow(t, 3);
    const double w19 = ik*w18;
    const double w20 = std::pow(u, 3);
    const double w21 = -w20;
    const double w22 = t*w3;
    const double w23 = t*w4;
    const double w24 = 15*w22 + 10*w23;
    const double w27 = -5*w6 - 5*w8;
    const double w28 = 2*w3;
    const double w29 = w0 + w28;
    const double w30 = 2*w15 + w4;
    const double w31 = 24*zt;
    const double w32 = w25*w31;
    const double w34 = 3*w15 + w27;
    const double w37 = w0*zt;
    const double w38 = u*w0;
    const double w39 = 10*w38;
    const double w40 = u*w15;
    const double w41 = 15*w40;
    const double w42 = 6*w38;
    const double w46 = 3*w3;
    const double w47 = w27 + w46;
    const double w49 = (-u + w7)/s;
    const double w50 = std::pow(t, 4);
    const double w51 = std::pow(u, 4);
    const double w52 = z*zt;
    const double w53 = 5*t*w20;
    const double w54 = u*w18;
    const double w55 = w0*w3;
    const double w56 = std::pow(u, 5);
    const double w58 = t*w51;
    const double w61 = w0*w20;
    const double w62 = w18*w3;
    o_mu_mu_mu_mu_mu[i] = w13*w2;
    o_mu_mu_mu_mu_sigma[i] = -w11*w17*w19;
    o_mu_mu_mu_mu_nu[i] = w2*w26*(20*ik*u*w0 + 15*u*w0 - w18 - 2*w19 - w21 - w24);
    o_mu_mu_mu_sigma_sigma[i] = -w10*w12*w15;
    o_mu_mu_mu_sigma_nu[i] = w2*w32*(-w27 - w29 - w30);
    o_mu_mu_mu_nu_nu[i] = w33*w36*w37*(-w29 - w34);
    o_mu_mu_sigma_sigma_sigma[i] = w11*w17*w18;
    o_mu_mu_sigma_sigma_nu[i] = w15*w26*(15*ik*t*w3 + ik*w18 - ik*w20 + 20*t*w3 - 2*w20 - w39 - w41);
    o_mu_mu_sigma_nu_nu[i] = w2*w43*(w19 - 3*w20 + w24 - 13*w40 - w42);
    o_mu_mu_nu_nu_nu[i] = 6*w18*w45*(w16 - 3*w6 - 7*w8);
    o_mu_sigma_sigma_sigma_sigma[i] = w0*w13;
    o_mu_sigma_sigma_sigma_nu[i] = w15*w32*(w15 + 2*w4 + w47);
    o_mu_sigma_sigma_nu_nu[i] = t*w1*w43*zt*(3*w19 + w21 + 13*w22 + 6*w23 - w39 - w41);
    o_mu_sigma_nu_nu_nu[i] = 12*w37*w45*(w15 + w28 - 2*w6 - 3*w8);
    o_mu_nu_nu_nu_nu[i] = w18*w31*w48*w49;
    o_sigma_sigma_sigma_sigma_sigma[i] = w11*(-std::pow(t, 5) + 10*u*w50 - w14*w18 + 15*w50*w52 + w51*w52 + w52*w53 + 5*w52*w54 + 10*w52*w55);
    o_sigma_sigma_sigma_sigma_nu[i] = u*w26*(2*ik*w0*w3 + 10*ik*w50 - 20*u*w19 - w51 - w53 - 35*w54 - 7*w55);
    o_sigma_sigma_sigma_nu_nu[i] = w36*w40*(w30 + w47);
    o_sigma_sigma_nu_nu_nu[i] = w1*w42*w44*(w34 + w46);
    o_sigma_nu_nu_nu_nu[i] = 24*w33*w49*w54;
    o_nu_nu_nu_nu_nu[i] = V7 + u*w50*w59 - 3*w48*w56 - 15*w48*w58 + 12*w56*w57 + w58*w59 - w60*w61 - w60*w62 + w61*w63 + w62*w63;
  });
  return List::create(Named("mu_mu_mu_mu_mu") = o_mu_mu_mu_mu_mu, Named("mu_mu_mu_mu_sigma") = o_mu_mu_mu_mu_sigma, Named("mu_mu_mu_mu_nu") = o_mu_mu_mu_mu_nu, Named("mu_mu_mu_sigma_sigma") = o_mu_mu_mu_sigma_sigma, Named("mu_mu_mu_sigma_nu") = o_mu_mu_mu_sigma_nu, Named("mu_mu_mu_nu_nu") = o_mu_mu_mu_nu_nu, Named("mu_mu_sigma_sigma_sigma") = o_mu_mu_sigma_sigma_sigma, Named("mu_mu_sigma_sigma_nu") = o_mu_mu_sigma_sigma_nu, Named("mu_mu_sigma_nu_nu") = o_mu_mu_sigma_nu_nu, Named("mu_mu_nu_nu_nu") = o_mu_mu_nu_nu_nu, Named("mu_sigma_sigma_sigma_sigma") = o_mu_sigma_sigma_sigma_sigma, Named("mu_sigma_sigma_sigma_nu") = o_mu_sigma_sigma_sigma_nu, Named("mu_sigma_sigma_nu_nu") = o_mu_sigma_sigma_nu_nu, Named("mu_sigma_nu_nu_nu") = o_mu_sigma_nu_nu_nu, Named("mu_nu_nu_nu_nu") = o_mu_nu_nu_nu_nu, Named("sigma_sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma_sigma, Named("sigma_sigma_sigma_sigma_nu") = o_sigma_sigma_sigma_sigma_nu, Named("sigma_sigma_sigma_nu_nu") = o_sigma_sigma_sigma_nu_nu, Named("sigma_sigma_nu_nu_nu") = o_sigma_sigma_nu_nu_nu, Named("sigma_nu_nu_nu_nu") = o_sigma_nu_nu_nu_nu, Named("nu_nu_nu_nu_nu") = o_nu_nu_nu_nu_nu);
}

// first derivatives of the expected Hessian
// [[Rcpp::export]]
List student_t1_dexpected1_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu_mu(n);
  NumericVector o_mu_mu_sigma(n);
  NumericVector o_mu_mu_nu(n);
  NumericVector o_sigma_sigma_mu(n);
  NumericVector o_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_nu(n);
  NumericVector o_nu_nu_mu(n);
  NumericVector o_nu_nu_sigma(n);
  NumericVector o_nu_nu_nu(n);
  NumericVector o_mu_sigma_mu(n);
  NumericVector o_mu_sigma_sigma(n);
  NumericVector o_mu_sigma_nu(n);
  NumericVector o_mu_nu_mu(n);
  NumericVector o_mu_nu_sigma(n);
  NumericVector o_mu_nu_nu(n);
  NumericVector o_sigma_nu_mu(n);
  NumericVector o_sigma_nu_sigma(n);
  NumericVector o_sigma_nu_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, V8, w0, w1, w2, w3, w4; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double V8 = t1_V8(v);
    const double w0 = 1/(std::pow(s, 3)*(3*ik + 1));
    const double w1 = std::pow(ik, 2);
    const double w2 = std::pow(s, -2);
    const double w3 = w1*w2/(6*ik + 9*w1 + 1);
    const double w4 = std::pow(ik, 3);
    return Prm{m, s, ik, V8, w0, w1, w2, w3, w4};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double V8 = P.V8; (void) V8;
    const double w0 = P.w0; (void) w0;
    const double w1 = P.w1; (void) w1;
    const double w2 = P.w2; (void) w2;
    const double w3 = P.w3; (void) w3;
    const double w4 = P.w4; (void) w4;
    o_mu_mu_mu[i] = 0;
    o_mu_mu_sigma[i] = 2*w0*(ik + 1);
    o_mu_mu_nu[i] = -2*w3;
    o_sigma_sigma_mu[i] = 0;
    o_sigma_sigma_sigma[i] = 4*w0;
    o_sigma_sigma_nu[i] = -6*w3;
    o_nu_nu_mu[i] = 0;
    o_nu_nu_sigma[i] = 0;
    o_nu_nu_nu[i] = V8;
    o_mu_sigma_mu[i] = 0;
    o_mu_sigma_sigma[i] = 0;
    o_mu_sigma_nu[i] = 0;
    o_mu_nu_mu[i] = 0;
    o_mu_nu_sigma[i] = 0;
    o_mu_nu_nu[i] = 0;
    o_sigma_nu_mu[i] = 0;
    o_sigma_nu_sigma[i] = -2*w1*w2/(4*ik + 3*w1 + 1);
    o_sigma_nu_nu[i] = -4*w4*(2*ik + 1)/(s*(9*std::pow(ik, 4) + 8*ik + 22*w1 + 24*w4 + 1));
  });
  return List::create(Named("mu_mu_mu") = o_mu_mu_mu, Named("mu_mu_sigma") = o_mu_mu_sigma, Named("mu_mu_nu") = o_mu_mu_nu, Named("sigma_sigma_mu") = o_sigma_sigma_mu, Named("sigma_sigma_sigma") = o_sigma_sigma_sigma, Named("sigma_sigma_nu") = o_sigma_sigma_nu, Named("nu_nu_mu") = o_nu_nu_mu, Named("nu_nu_sigma") = o_nu_nu_sigma, Named("nu_nu_nu") = o_nu_nu_nu, Named("mu_sigma_mu") = o_mu_sigma_mu, Named("mu_sigma_sigma") = o_mu_sigma_sigma, Named("mu_sigma_nu") = o_mu_sigma_nu, Named("mu_nu_mu") = o_mu_nu_mu, Named("mu_nu_sigma") = o_mu_nu_sigma, Named("mu_nu_nu") = o_mu_nu_nu, Named("sigma_nu_mu") = o_sigma_nu_mu, Named("sigma_nu_sigma") = o_sigma_nu_sigma, Named("sigma_nu_nu") = o_sigma_nu_nu);
}

// second derivatives of the expected Hessian
// [[Rcpp::export]]
List student_t1_dexpected2_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu_mu_mu(n);
  NumericVector o_mu_mu_sigma_sigma(n);
  NumericVector o_mu_mu_nu_nu(n);
  NumericVector o_mu_mu_mu_sigma(n);
  NumericVector o_mu_mu_mu_nu(n);
  NumericVector o_mu_mu_sigma_nu(n);
  NumericVector o_sigma_sigma_mu_mu(n);
  NumericVector o_sigma_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_nu_nu(n);
  NumericVector o_sigma_sigma_mu_sigma(n);
  NumericVector o_sigma_sigma_mu_nu(n);
  NumericVector o_sigma_sigma_sigma_nu(n);
  NumericVector o_nu_nu_mu_mu(n);
  NumericVector o_nu_nu_sigma_sigma(n);
  NumericVector o_nu_nu_nu_nu(n);
  NumericVector o_nu_nu_mu_sigma(n);
  NumericVector o_nu_nu_mu_nu(n);
  NumericVector o_nu_nu_sigma_nu(n);
  NumericVector o_mu_sigma_mu_mu(n);
  NumericVector o_mu_sigma_sigma_sigma(n);
  NumericVector o_mu_sigma_nu_nu(n);
  NumericVector o_mu_sigma_mu_sigma(n);
  NumericVector o_mu_sigma_mu_nu(n);
  NumericVector o_mu_sigma_sigma_nu(n);
  NumericVector o_mu_nu_mu_mu(n);
  NumericVector o_mu_nu_sigma_sigma(n);
  NumericVector o_mu_nu_nu_nu(n);
  NumericVector o_mu_nu_mu_sigma(n);
  NumericVector o_mu_nu_mu_nu(n);
  NumericVector o_mu_nu_sigma_nu(n);
  NumericVector o_sigma_nu_mu_mu(n);
  NumericVector o_sigma_nu_sigma_sigma(n);
  NumericVector o_sigma_nu_nu_nu(n);
  NumericVector o_sigma_nu_mu_sigma(n);
  NumericVector o_sigma_nu_mu_nu(n);
  NumericVector o_sigma_nu_sigma_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, V9, w0, w1, w2, w3, w4, w5, w6, w7, w8; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double V9 = t1_V9(v);
    const double w0 = 1/(std::pow(s, 4)*(3*ik + 1));
    const double w1 = std::pow(ik, 3);
    const double w2 = std::pow(s, -2);
    const double w3 = std::pow(ik, 2);
    const double w4 = w1*w2/(9*ik + 27*w1 + 27*w3 + 1);
    const double w5 = std::pow(s, -3);
    const double w6 = w3*w5/(6*ik + 9*w3 + 1);
    const double w7 = std::pow(ik, 4);
    const double w8 = 12*ik;
    return Prm{m, s, ik, V9, w0, w1, w2, w3, w4, w5, w6, w7, w8};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double V9 = P.V9; (void) V9;
    const double w0 = P.w0; (void) w0;
    const double w1 = P.w1; (void) w1;
    const double w2 = P.w2; (void) w2;
    const double w3 = P.w3; (void) w3;
    const double w4 = P.w4; (void) w4;
    const double w5 = P.w5; (void) w5;
    const double w6 = P.w6; (void) w6;
    const double w7 = P.w7; (void) w7;
    const double w8 = P.w8; (void) w8;
    o_mu_mu_mu_mu[i] = 0;
    o_mu_mu_sigma_sigma[i] = -6*w0*(ik + 1);
    o_mu_mu_nu_nu[i] = 4*w4;
    o_mu_mu_mu_sigma[i] = 0;
    o_mu_mu_mu_nu[i] = 0;
    o_mu_mu_sigma_nu[i] = 4*w6;
    o_sigma_sigma_mu_mu[i] = 0;
    o_sigma_sigma_sigma_sigma[i] = -12*w0;
    o_sigma_sigma_nu_nu[i] = 12*w4;
    o_sigma_sigma_mu_sigma[i] = 0;
    o_sigma_sigma_mu_nu[i] = 0;
    o_sigma_sigma_sigma_nu[i] = 12*w6;
    o_nu_nu_mu_mu[i] = 0;
    o_nu_nu_sigma_sigma[i] = 0;
    o_nu_nu_nu_nu[i] = V9;
    o_nu_nu_mu_sigma[i] = 0;
    o_nu_nu_mu_nu[i] = 0;
    o_nu_nu_sigma_nu[i] = 0;
    o_mu_sigma_mu_mu[i] = 0;
    o_mu_sigma_sigma_sigma[i] = 0;
    o_mu_sigma_nu_nu[i] = 0;
    o_mu_sigma_mu_sigma[i] = 0;
    o_mu_sigma_mu_nu[i] = 0;
    o_mu_sigma_sigma_nu[i] = 0;
    o_mu_nu_mu_mu[i] = 0;
    o_mu_nu_sigma_sigma[i] = 0;
    o_mu_nu_nu_nu[i] = 0;
    o_mu_nu_mu_sigma[i] = 0;
    o_mu_nu_mu_nu[i] = 0;
    o_mu_nu_sigma_nu[i] = 0;
    o_sigma_nu_mu_mu[i] = 0;
    o_sigma_nu_sigma_sigma[i] = 4*w3*w5/(4*ik + 3*w3 + 1);
    o_sigma_nu_nu_nu[i] = 4*w7*(13*w3 + w8 + 3)/(s*(27*std::pow(ik, 6) + 108*std::pow(ik, 5) + 136*w1 + 57*w3 + 171*w7 + w8 + 1));
    o_sigma_nu_mu_sigma[i] = 0;
    o_sigma_nu_mu_nu[i] = 0;
    o_sigma_nu_sigma_nu[i] = 4*w1*w2*(2*ik + 1)/(8*ik + 24*w1 + 22*w3 + 9*w7 + 1);
  });
  return List::create(Named("mu_mu_mu_mu") = o_mu_mu_mu_mu, Named("mu_mu_sigma_sigma") = o_mu_mu_sigma_sigma, Named("mu_mu_nu_nu") = o_mu_mu_nu_nu, Named("mu_mu_mu_sigma") = o_mu_mu_mu_sigma, Named("mu_mu_mu_nu") = o_mu_mu_mu_nu, Named("mu_mu_sigma_nu") = o_mu_mu_sigma_nu, Named("sigma_sigma_mu_mu") = o_sigma_sigma_mu_mu, Named("sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma, Named("sigma_sigma_nu_nu") = o_sigma_sigma_nu_nu, Named("sigma_sigma_mu_sigma") = o_sigma_sigma_mu_sigma, Named("sigma_sigma_mu_nu") = o_sigma_sigma_mu_nu, Named("sigma_sigma_sigma_nu") = o_sigma_sigma_sigma_nu, Named("nu_nu_mu_mu") = o_nu_nu_mu_mu, Named("nu_nu_sigma_sigma") = o_nu_nu_sigma_sigma, Named("nu_nu_nu_nu") = o_nu_nu_nu_nu, Named("nu_nu_mu_sigma") = o_nu_nu_mu_sigma, Named("nu_nu_mu_nu") = o_nu_nu_mu_nu, Named("nu_nu_sigma_nu") = o_nu_nu_sigma_nu, Named("mu_sigma_mu_mu") = o_mu_sigma_mu_mu, Named("mu_sigma_sigma_sigma") = o_mu_sigma_sigma_sigma, Named("mu_sigma_nu_nu") = o_mu_sigma_nu_nu, Named("mu_sigma_mu_sigma") = o_mu_sigma_mu_sigma, Named("mu_sigma_mu_nu") = o_mu_sigma_mu_nu, Named("mu_sigma_sigma_nu") = o_mu_sigma_sigma_nu, Named("mu_nu_mu_mu") = o_mu_nu_mu_mu, Named("mu_nu_sigma_sigma") = o_mu_nu_sigma_sigma, Named("mu_nu_nu_nu") = o_mu_nu_nu_nu, Named("mu_nu_mu_sigma") = o_mu_nu_mu_sigma, Named("mu_nu_mu_nu") = o_mu_nu_mu_nu, Named("mu_nu_sigma_nu") = o_mu_nu_sigma_nu, Named("sigma_nu_mu_mu") = o_sigma_nu_mu_mu, Named("sigma_nu_sigma_sigma") = o_sigma_nu_sigma_sigma, Named("sigma_nu_nu_nu") = o_sigma_nu_nu_nu, Named("sigma_nu_mu_sigma") = o_sigma_nu_mu_sigma, Named("sigma_nu_mu_nu") = o_sigma_nu_mu_nu, Named("sigma_nu_sigma_nu") = o_sigma_nu_sigma_nu);
}

// 1 in y and 1 in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t1_cross_y_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu(n);
  NumericVector o_sigma(n);
  NumericVector o_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    return Prm{m, s, ik};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    // u = q/(1 + q) in [0, 1], and zt = z t, where q itself may overflow
    const double u = (q < 1.0) ? q * t : 1.0 / (1.0 + 1.0 / q);
    const double zt = (std::fabs(z) <= 1.0) ? z * t : 1.0 / (1.0 / z + z * ik);
    (void) t; (void) q; (void) u; (void) zt;
    const double w0 = ik*t - u;
    const double w1 = t/std::pow(s, 2);
    o_mu[i] = w1*(-ik*u + t + w0);
    o_sigma[i] = 2*w1*zt*(ik + 1);
    o_nu[i] = ik*w0*zt/s;
  });
  return List::create(Named("mu") = o_mu, Named("sigma") = o_sigma, Named("nu") = o_nu);
}

// 2 in y and 1 in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t1_cross2_y_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu(n);
  NumericVector o_sigma(n);
  NumericVector o_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, w2; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double w2 = 2/std::pow(s, 3);
    return Prm{m, s, ik, w2};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double w2 = P.w2; (void) w2;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    // u = q/(1 + q) in [0, 1], and zt = z t, where q itself may overflow
    const double u = (q < 1.0) ? q * t : 1.0 / (1.0 + 1.0 / q);
    const double zt = (std::fabs(z) <= 1.0) ? z * t : 1.0 / (1.0 / z + z * ik);
    (void) t; (void) q; (void) u; (void) zt;
    const double w0 = 3*t;
    const double w1 = ik*t;
    const double w3 = std::pow(t, 2);
    const double w4 = ik*u;
    o_mu[i] = w1*w2*zt*(ik*u - ik*w0 + u - w0);
    o_sigma[i] = w2*w3*(t - 3*u + w1 - 3*w4);
    o_nu[i] = w1*(ik*w3 + std::pow(u, 2) - u*w0 - w0*w4)/std::pow(s, 2);
  });
  return List::create(Named("mu") = o_mu, Named("sigma") = o_sigma, Named("nu") = o_nu);
}

// 1 in y and 2 in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t1_grad_y_hess_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu(n);
  NumericVector o_sigma_sigma(n);
  NumericVector o_nu_nu(n);
  NumericVector o_mu_sigma(n);
  NumericVector o_mu_nu(n);
  NumericVector o_sigma_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, w5; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double w5 = 2/std::pow(s, 3);
    return Prm{m, s, ik, w5};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double w5 = P.w5; (void) w5;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    // u = q/(1 + q) in [0, 1], and zt = z t, where q itself may overflow
    const double u = (q < 1.0) ? q * t : 1.0 / (1.0 + 1.0 / q);
    const double zt = (std::fabs(z) <= 1.0) ? z * t : 1.0 / (1.0 / z + z * ik);
    (void) t; (void) q; (void) u; (void) zt;
    const double w0 = 3*t;
    const double w1 = -u;
    const double w2 = ik*u;
    const double w3 = ik*w0 + w0 + w1 - w2;
    const double w4 = ik*t;
    const double w6 = w5*zt;
    const double w7 = 2*zt;
    const double w8 = std::pow(t, 2);
    const double w9 = w4/std::pow(s, 2);
    o_mu_mu[i] = w3*w4*w6;
    o_sigma_sigma[i] = -t*w3*w6;
    o_nu_nu[i] = std::pow(ik, 2)*t*w7*(-w1 - w4)/s;
    o_mu_sigma[i] = w5*w8*(3*ik*u - t + 3*u - w4);
    o_mu_nu[i] = w9*(3*ik*t*u - ik*w8 + 3*t*u - std::pow(u, 2));
    o_sigma_nu[i] = w7*w9*(2*u + w2 - w4);
  });
  return List::create(Named("mu_mu") = o_mu_mu, Named("sigma_sigma") = o_sigma_sigma, Named("nu_nu") = o_nu_nu, Named("mu_sigma") = o_mu_sigma, Named("mu_nu") = o_mu_nu, Named("sigma_nu") = o_sigma_nu);
}

// 2 in y and 2 in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t1_hess_y_hess_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu(n);
  NumericVector o_sigma_sigma(n);
  NumericVector o_nu_nu(n);
  NumericVector o_mu_sigma(n);
  NumericVector o_mu_nu(n);
  NumericVector o_sigma_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, w0, w11, w12; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double w0 = std::pow(s, -4);
    const double w11 = 2*std::pow(ik, 2);
    const double w12 = std::pow(s, -3);
    return Prm{m, s, ik, w0, w11, w12};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double w0 = P.w0; (void) w0;
    const double w11 = P.w11; (void) w11;
    const double w12 = P.w12; (void) w12;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    // u = q/(1 + q) in [0, 1], and zt = z t, where q itself may overflow
    const double u = (q < 1.0) ? q * t : 1.0 / (1.0 + 1.0 / q);
    const double zt = (std::fabs(z) <= 1.0) ? z * t : 1.0 / (1.0 / z + z * ik);
    (void) t; (void) q; (void) u; (void) zt;
    const double w1 = std::pow(t, 2);
    const double w2 = std::pow(u, 2);
    const double w3 = ik*w2;
    const double w4 = ik*t;
    const double w5 = -6*u*w4 + w2;
    const double w6 = ik*w1;
    const double w7 = t*u;
    const double w8 = w6 - 6*w7;
    const double w9 = w1 + w3 + w5 + w8;
    const double w10 = 6*w6;
    o_mu_mu[i] = w0*w10*w9;
    o_sigma_sigma[i] = -6*w0*w1*w9;
    o_nu_nu[i] = w1*w11*(5*ik*t*u + 3*t*u - 3*w2 - w6)/std::pow(s, 2);
    o_mu_sigma[i] = 24*w0*w6*zt*(-ik*u + t - u + w4);
    o_mu_nu[i] = t*w11*w12*zt*(3*w1 + w10 + w5 - 8*w7);
    o_sigma_nu[i] = 2*w12*w6*(8*ik*t*u - 6*w2 - 3*w3 - w8);
  });
  return List::create(Named("mu_mu") = o_mu_mu, Named("sigma_sigma") = o_sigma_sigma, Named("nu_nu") = o_nu_nu, Named("mu_sigma") = o_mu_sigma, Named("mu_nu") = o_mu_nu, Named("sigma_nu") = o_sigma_nu);
}

// order 1 of log f in y
// [[Rcpp::export]]
List student_t1_dy1_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_y(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    return Prm{m, s, ik};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    // u = q/(1 + q) in [0, 1], and zt = z t, where q itself may overflow
    const double u = (q < 1.0) ? q * t : 1.0 / (1.0 + 1.0 / q);
    const double zt = (std::fabs(z) <= 1.0) ? z * t : 1.0 / (1.0 / z + z * ik);
    (void) t; (void) q; (void) u; (void) zt;
    o_y[i] = -zt*(ik + 1)/s;
  });
  return List::create(Named("y") = o_y);
}

// order 2 of log f in y
// [[Rcpp::export]]
List student_t1_dy2_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_y(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    return Prm{m, s, ik};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    // u = q/(1 + q) in [0, 1], and zt = z t, where q itself may overflow
    const double u = (q < 1.0) ? q * t : 1.0 / (1.0 + 1.0 / q);
    const double zt = (std::fabs(z) <= 1.0) ? z * t : 1.0 / (1.0 / z + z * ik);
    (void) t; (void) q; (void) u; (void) zt;
    o_y[i] = t*(-ik*t + ik*u - t + u)/std::pow(s, 2);
  });
  return List::create(Named("y") = o_y);
}

// order 3 of log f in y
// [[Rcpp::export]]
List student_t1_dy3_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_y(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    return Prm{m, s, ik};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    // u = q/(1 + q) in [0, 1], and zt = z t, where q itself may overflow
    const double u = (q < 1.0) ? q * t : 1.0 / (1.0 + 1.0 / q);
    const double zt = (std::fabs(z) <= 1.0) ? z * t : 1.0 / (1.0 / z + z * ik);
    (void) t; (void) q; (void) u; (void) zt;
    const double w0 = 3*t;
    o_y[i] = 2*ik*t*zt*(-ik*u + ik*w0 - u + w0)/std::pow(s, 3);
  });
  return List::create(Named("y") = o_y);
}

// order 4 of log f in y
// [[Rcpp::export]]
List student_t1_dy4_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_y(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    return Prm{m, s, ik};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    // u = q/(1 + q) in [0, 1], and zt = z t, where q itself may overflow
    const double u = (q < 1.0) ? q * t : 1.0 / (1.0 + 1.0 / q);
    const double zt = (std::fabs(z) <= 1.0) ? z * t : 1.0 / (1.0 / z + z * ik);
    (void) t; (void) q; (void) u; (void) zt;
    const double w0 = std::pow(t, 2);
    const double w1 = std::pow(u, 2);
    const double w2 = 6*t*u;
    const double w3 = ik*w0;
    o_y[i] = 6*w3*(ik*w1 - ik*w2 + w0 + w1 - w2 + w3)/std::pow(s, 4);
  });
  return List::create(Named("y") = o_y);
}
