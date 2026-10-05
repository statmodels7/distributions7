#include <Rcpp.h>
#include <cmath>
#include "d7_par.h"
#include "pt_student_t2.h"
using namespace Rcpp;

// The Student t in its location mu, standard deviation s and degrees of
// freedom nu > 2. With z = (y - mu)/s, k = nu - 2, q = z^2/k and
// t = 1/(1 + q), the log-density is
//   c(nu) - log s - log(pi)/2 - (k + 3)/2 log(1 + q),
//   c(nu) = lgamma((nu+1)/2) - lgamma(nu/2) - log(nu - 2)/2.
// Every component below is a closed form derived offline with sympy
// (stabilita/gen_student_t2.py, T_FAMILY=t2), one exported function per
// order and per surface, and is checked numerically in the tests.
//
// Every derivative in nu vanishes as nu grows and is a difference of terms
// agreeing to leading order. The data part is cancelled symbolically and
// written in q, z, 1/k and t, which neither cancel nor overflow at the nu
// the link can produce; the score's logarithm enters as
// D(q) = q/(1+q) - log1p(q). Each quantity of nu alone that carries a
// polygamma, t2_V*, has two branches: the direct form below nu = 20, and
// above it the asymptotic series in h = 1/nu, exact to h^40, from
// Stirling's series. The expectations are closed, t being Beta(nu/2, 1/2)
// under the model.

// D(q) and the quantities of nu alone, t2_*, are in pt_student_t2.h.

// order 1 of log f in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t2_gradient_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
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
    const double ik = 1.0 / (v - 2.0);
    const double V0 = t2_V0(v);
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
    const double DQ = d7::student_t2_DQ(z, ik);
    o_mu[i] = d7::student_t2_score_mu(z, s, ik);
    o_sigma[i] = d7::student_t2_score_sigma(z, s, ik);
    o_nu[i] = d7::student_t2_score_nu(z, ik, V0, DQ);
  });
  return List::create(Named("mu") = o_mu, Named("sigma") = o_sigma, Named("nu") = o_nu);
}

// order 2 of log f in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t2_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu(n);
  NumericVector o_sigma_sigma(n);
  NumericVector o_nu_nu(n);
  NumericVector o_mu_sigma(n);
  NumericVector o_mu_nu(n);
  NumericVector o_sigma_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, V1, w0, w2, w7; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 2.0);
    const double V1 = t2_V1(v);
    const double w0 = 3*ik;
    const double w2 = std::pow(s, -2);
    const double w7 = std::pow(ik, 2);
    return Prm{m, s, ik, V1, w0, w2, w7};
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
    const double w0 = P.w0; (void) w0;
    const double w2 = P.w2; (void) w2;
    const double w7 = P.w7; (void) w7;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    // u = q/(1 + q) in [0, 1], and zt = z t, where q itself may overflow
    const double u = (q < 1.0) ? q * t : 1.0 / (1.0 + 1.0 / q);
    const double zt = (std::fabs(z) <= 1.0) ? z * t : 1.0 / (1.0 / z + z * ik);
    (void) t; (void) q; (void) u; (void) zt;
    const double w1 = t*w0 - u;
    const double w3 = t*w2;
    const double w8 = -w1/s;
    o_mu_mu[i] = d7::student_t2_hess_mu_mu(z, s, ik);
    o_sigma_sigma[i] = d7::student_t2_hess_sigma_sigma(z, s, ik);
    o_nu_nu[i] = d7::student_t2_hess_nu_nu(z, ik, V1);
    o_mu_sigma[i] = -2*w3*zt*(w0 + 1);
    o_mu_nu[i] = ik*w8*zt;
    o_sigma_nu[i] = u*w8;
  });
  return List::create(Named("mu_mu") = o_mu_mu, Named("sigma_sigma") = o_sigma_sigma, Named("nu_nu") = o_nu_nu, Named("mu_sigma") = o_mu_sigma, Named("mu_nu") = o_mu_nu, Named("sigma_nu") = o_sigma_nu);
}

// order 2, expected
// [[Rcpp::export]]
List student_t2_expected_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu(n);
  NumericVector o_sigma_sigma(n);
  NumericVector o_nu_nu(n);
  NumericVector o_mu_sigma(n);
  NumericVector o_mu_nu(n);
  NumericVector o_sigma_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, V2, w0, w1, w2; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 2.0);
    const double V2 = t2_V2(v);
    const double w0 = std::pow(ik, 2);
    const double w1 = 5*ik + 1;
    const double w2 = 1/(std::pow(s, 2)*w1);
    return Prm{m, s, ik, V2, w0, w1, w2};
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
    const double w2 = P.w2; (void) w2;
    o_mu_mu[i] = d7::student_t2_expected_mu_mu(s, ik);
    o_sigma_sigma[i] = d7::student_t2_expected_sigma_sigma(s, ik);
    o_nu_nu[i] = d7::student_t2_expected_nu_nu(V2);
    o_mu_sigma[i] = 0;
    o_mu_nu[i] = 0;
    o_sigma_nu[i] = -6*std::pow(ik, 3)/(s*(8*ik + 15*w0 + 1));
  });
  return List::create(Named("mu_mu") = o_mu_mu, Named("sigma_sigma") = o_sigma_sigma, Named("nu_nu") = o_nu_nu, Named("mu_sigma") = o_mu_sigma, Named("mu_nu") = o_mu_nu, Named("sigma_nu") = o_sigma_nu);
}

// order 3 of log f in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t2_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
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
  struct Prm { double m, s, ik, V3, w2, w8, w14, w18, w24, w25; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 2.0);
    const double V3 = t2_V3(v);
    const double w2 = 9*ik;
    const double w8 = 2/std::pow(s, 3);
    const double w14 = std::pow(s, -2);
    const double w18 = std::pow(ik, 2);
    const double w24 = std::pow(ik, 3);
    const double w25 = 9*w24;
    return Prm{m, s, ik, V3, w2, w8, w14, w18, w24, w25};
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
    const double w2 = P.w2; (void) w2;
    const double w8 = P.w8; (void) w8;
    const double w14 = P.w14; (void) w14;
    const double w18 = P.w18; (void) w18;
    const double w24 = P.w24; (void) w24;
    const double w25 = P.w25; (void) w25;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    // u = q/(1 + q) in [0, 1], and zt = z t, where q itself may overflow
    const double u = (q < 1.0) ? q * t : 1.0 / (1.0 + 1.0 / q);
    const double zt = (std::fabs(z) <= 1.0) ? z * t : 1.0 / (1.0 / z + z * ik);
    (void) t; (void) q; (void) u; (void) zt;
    const double w0 = 3*t;
    const double w1 = -u;
    const double w3 = t*w2;
    const double w4 = 3*u;
    const double w5 = ik*w4;
    const double w6 = w0 + w1 + w3 - w5;
    const double w7 = ik*t;
    const double w9 = w8*zt;
    const double w10 = std::pow(t, 2);
    const double w11 = ik*w0;
    const double w12 = std::pow(u, 2);
    const double w13 = u*w0;
    const double w15 = w14*w7;
    const double w16 = 2*u;
    const double w17 = 2*zt;
    const double w19 = (w1 + w11)/s;
    const double w20 = std::pow(u, 3);
    const double w21 = t*w12;
    const double w22 = u*w10;
    const double w23 = z*zt;
    o_mu_mu_mu[i] = -w6*w7*w9;
    o_mu_mu_sigma[i] = w10*w8*(t - u*w2 + w11 - w4);
    o_mu_mu_nu[i] = w15*(3*ik*w10 - u*w3 + w12 - w13);
    o_mu_sigma_sigma[i] = t*w6*w9;
    o_mu_sigma_nu[i] = w15*w17*(3*ik*t - w16 - w5);
    o_mu_nu_nu[i] = t*w17*w18*w19;
    o_sigma_sigma_sigma[i] = w8*(-std::pow(t, 3) + 6*w10*w23 + w12*w23 + w13*w23 + 2*w20 + 6*w21 + 15*w22);
    o_sigma_sigma_nu[i] = u*w14*(9*ik*w10 - 5*t*u - u*w11 - w12);
    o_sigma_nu_nu[i] = w16*w19*w7;
    o_nu_nu_nu[i] = V3 - 1.0/2.0*w18*w20 - 3.0/2.0*w18*w21 + 3*w20*w24 + w21*w25 + w22*w25;
  });
  return List::create(Named("mu_mu_mu") = o_mu_mu_mu, Named("mu_mu_sigma") = o_mu_mu_sigma, Named("mu_mu_nu") = o_mu_mu_nu, Named("mu_sigma_sigma") = o_mu_sigma_sigma, Named("mu_sigma_nu") = o_mu_sigma_nu, Named("mu_nu_nu") = o_mu_nu_nu, Named("sigma_sigma_sigma") = o_sigma_sigma_sigma, Named("sigma_sigma_nu") = o_sigma_sigma_nu, Named("sigma_nu_nu") = o_sigma_nu_nu, Named("nu_nu_nu") = o_nu_nu_nu);
}

// order 3, expected
// [[Rcpp::export]]
List student_t2_deriv3_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
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
  struct Prm { double m, s, ik, V4, w0, w1, w2, w3, w4, w5; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 2.0);
    const double V4 = t2_V4(v);
    const double w0 = std::pow(ik, 2);
    const double w1 = std::pow(ik, 3);
    const double w2 = 1.0/(12*ik + 35*w0 + 1);
    const double w3 = 2*w2/std::pow(s, 3);
    const double w4 = 2*ik + 1;
    const double w5 = w2/std::pow(s, 2);
    return Prm{m, s, ik, V4, w0, w1, w2, w3, w4, w5};
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
    o_mu_mu_mu[i] = 0;
    o_mu_mu_sigma[i] = w3*(9*ik + 26*w0 + 24*w1 + 1);
    o_mu_mu_nu[i] = 12*w1*w4*w5;
    o_mu_sigma_sigma[i] = 0;
    o_mu_sigma_nu[i] = 0;
    o_mu_nu_nu[i] = 0;
    o_sigma_sigma_sigma[i] = w3*(33*ik + 46*w0 + 5);
    o_sigma_sigma_nu[i] = 6*w0*w5*(3*ik - 1);
    o_sigma_nu_nu[i] = 24*std::pow(ik, 4)*w4/(s*(15*ik + 71*w0 + 105*w1 + 1));
    o_nu_nu_nu[i] = V4;
  });
  return List::create(Named("mu_mu_mu") = o_mu_mu_mu, Named("mu_mu_sigma") = o_mu_mu_sigma, Named("mu_mu_nu") = o_mu_mu_nu, Named("mu_sigma_sigma") = o_mu_sigma_sigma, Named("mu_sigma_nu") = o_mu_sigma_nu, Named("mu_nu_nu") = o_mu_nu_nu, Named("sigma_sigma_sigma") = o_sigma_sigma_sigma, Named("sigma_sigma_nu") = o_sigma_sigma_nu, Named("sigma_nu_nu") = o_sigma_nu_nu, Named("nu_nu_nu") = o_nu_nu_nu);
}

// order 4 of log f in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t2_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
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
  struct Prm { double m, s, ik, V5, w1, w4, w11, w12, w15, w16, w18, w19, w23, w24, w28, w29, w40, w41; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 2.0);
    const double V5 = t2_V5(v);
    const double w1 = 18*ik;
    const double w4 = 3*ik;
    const double w11 = std::pow(s, -4);
    const double w12 = 6*w11;
    const double w15 = 24*w11;
    const double w16 = std::pow(ik, 2);
    const double w18 = std::pow(s, -3);
    const double w19 = 2*w18;
    const double w23 = std::pow(s, -2);
    const double w24 = w16*w23;
    const double w28 = std::pow(ik, 3);
    const double w29 = 6*w28;
    const double w40 = std::pow(ik, 4);
    const double w41 = 36*w40;
    return Prm{m, s, ik, V5, w1, w4, w11, w12, w15, w16, w18, w19, w23, w24, w28, w29, w40, w41};
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
    const double w1 = P.w1; (void) w1;
    const double w4 = P.w4; (void) w4;
    const double w11 = P.w11; (void) w11;
    const double w12 = P.w12; (void) w12;
    const double w15 = P.w15; (void) w15;
    const double w16 = P.w16; (void) w16;
    const double w18 = P.w18; (void) w18;
    const double w19 = P.w19; (void) w19;
    const double w23 = P.w23; (void) w23;
    const double w24 = P.w24; (void) w24;
    const double w28 = P.w28; (void) w28;
    const double w29 = P.w29; (void) w29;
    const double w40 = P.w40; (void) w40;
    const double w41 = P.w41; (void) w41;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    // u = q/(1 + q) in [0, 1], and zt = z t, where q itself may overflow
    const double u = (q < 1.0) ? q * t : 1.0 / (1.0 + 1.0 / q);
    const double zt = (std::fabs(z) <= 1.0) ? z * t : 1.0 / (1.0 / z + z * ik);
    (void) t; (void) q; (void) u; (void) zt;
    const double w0 = t*u;
    const double w2 = -w0*w1;
    const double w3 = std::pow(u, 2);
    const double w5 = w3*w4;
    const double w6 = std::pow(t, 2);
    const double w7 = 3*w6;
    const double w8 = ik*w7 + w3;
    const double w9 = -6*w0 + w2 + w5 + w6 + w8;
    const double w10 = ik*w6;
    const double w13 = t*w4 - u;
    const double w14 = t - u*w4 + w13;
    const double w17 = 18*w10;
    const double w20 = t*w19*zt;
    const double w21 = -2*t*u;
    const double w22 = 2*w3 + w5;
    const double w25 = w6*zt;
    const double w26 = 9*w10;
    const double w27 = 4*t;
    const double w30 = -w13/s;
    const double w31 = std::pow(u, 4);
    const double w32 = std::pow(u, 3);
    const double w33 = t*w32;
    const double w34 = std::pow(t, 3);
    const double w35 = u*w34;
    const double w36 = z*zt;
    const double w37 = w3*w6;
    const double w38 = w27*w3;
    const double w39 = u*w6;
    o_mu_mu_mu_mu[i] = w10*w12*w9;
    o_mu_mu_mu_sigma[i] = w10*w14*w15*zt;
    o_mu_mu_mu_nu[i] = w16*w20*(-8*w0 + w17 + w2 + w3 + w7);
    o_mu_mu_sigma_sigma[i] = -w12*w6*w9;
    o_mu_mu_sigma_nu[i] = 6*w10*w18*(8*ik*t*u - w10 - w21 - w22);
    o_mu_mu_nu_nu[i] = 6*w24*w6*(5*ik*t*u + t*u - w10 - w3);
    o_mu_sigma_sigma_sigma[i] = -w14*w15*w25;
    o_mu_sigma_sigma_nu[i] = ik*w20*(24*ik*t*u + 10*t*u - w22 - w26);
    o_mu_sigma_nu_nu[i] = w24*w27*zt*(6*ik*t*u - w21 - w8);
    o_mu_nu_nu_nu[i] = w25*w29*w30;
    o_sigma_sigma_sigma_sigma[i] = w12*(std::pow(t, 4) - 2*w31 - w32*w36 - 8*w33 - 10*w34*w36 - 26*w35 - w36*w38 - 5*w36*w39 - 9*w37);
    o_sigma_sigma_sigma_nu[i] = u*w19*(u*w17 - w1*w34 + w32 + w38 + 15*w39);
    o_sigma_sigma_nu_nu[i] = 2*ik*w0*w23*(9*ik*t*u + 5*t*u - w26 - w3);
    o_sigma_nu_nu_nu[i] = 6*w16*w30*w39;
    o_nu_nu_nu_nu[i] = V5 + w28*w31 + 4*w28*w33 + w29*w37 - 9*w31*w40 - w33*w41 - w35*w41 - 54*w37*w40;
  });
  return List::create(Named("mu_mu_mu_mu") = o_mu_mu_mu_mu, Named("mu_mu_mu_sigma") = o_mu_mu_mu_sigma, Named("mu_mu_mu_nu") = o_mu_mu_mu_nu, Named("mu_mu_sigma_sigma") = o_mu_mu_sigma_sigma, Named("mu_mu_sigma_nu") = o_mu_mu_sigma_nu, Named("mu_mu_nu_nu") = o_mu_mu_nu_nu, Named("mu_sigma_sigma_sigma") = o_mu_sigma_sigma_sigma, Named("mu_sigma_sigma_nu") = o_mu_sigma_sigma_nu, Named("mu_sigma_nu_nu") = o_mu_sigma_nu_nu, Named("mu_nu_nu_nu") = o_mu_nu_nu_nu, Named("sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma, Named("sigma_sigma_sigma_nu") = o_sigma_sigma_sigma_nu, Named("sigma_sigma_nu_nu") = o_sigma_sigma_nu_nu, Named("sigma_nu_nu_nu") = o_sigma_nu_nu_nu, Named("nu_nu_nu_nu") = o_nu_nu_nu_nu);
}

// order 4, expected
// [[Rcpp::export]]
List student_t2_deriv4_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
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
  struct Prm { double m, s, ik, V6, w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 2.0);
    const double V6 = t2_V6(v);
    const double w0 = 6*ik;
    const double w1 = std::pow(ik, 2);
    const double w2 = std::pow(ik, 3);
    const double w3 = 24*w2;
    const double w4 = 1/(std::pow(s, 4)*(16*ik + 63*w1 + 1));
    const double w5 = w4*(9*ik + 26*w1 + w3 + 1);
    const double w6 = 1.0/(21*ik + 143*w1 + 315*w2 + 1);
    const double w7 = 6*w1*w6/std::pow(s, 3);
    const double w8 = std::pow(ik, 4);
    const double w9 = w0 + 8*w1 + 1;
    const double w10 = w6/std::pow(s, 2);
    return Prm{m, s, ik, V6, w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10};
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
    const double w10 = P.w10; (void) w10;
    o_mu_mu_mu_mu[i] = w0*w5;
    o_mu_mu_mu_sigma[i] = 0;
    o_mu_mu_mu_nu[i] = 0;
    o_mu_mu_sigma_sigma[i] = -6*w5;
    o_mu_mu_sigma_nu[i] = -w7*(-3*ik + 10*w1 + w3 - 1);
    o_mu_mu_nu_nu[i] = -36*w10*w8*w9;
    o_mu_sigma_sigma_sigma[i] = 0;
    o_mu_sigma_sigma_nu[i] = 0;
    o_mu_sigma_nu_nu[i] = 0;
    o_mu_nu_nu_nu[i] = 0;
    o_sigma_sigma_sigma_sigma[i] = -18*w4*(23*ik + 34*w1 + 3);
    o_sigma_sigma_sigma_nu[i] = w7*(29*ik - 48*w1 + 9);
    o_sigma_sigma_nu_nu[i] = -12*w10*w2*(4*ik + 12*w1 - 1);
    o_sigma_nu_nu_nu[i] = -108*std::pow(ik, 5)*w9/(s*(24*ik + 206*w1 + 744*w2 + 945*w8 + 1));
    o_nu_nu_nu_nu[i] = V6;
  });
  return List::create(Named("mu_mu_mu_mu") = o_mu_mu_mu_mu, Named("mu_mu_mu_sigma") = o_mu_mu_mu_sigma, Named("mu_mu_mu_nu") = o_mu_mu_mu_nu, Named("mu_mu_sigma_sigma") = o_mu_mu_sigma_sigma, Named("mu_mu_sigma_nu") = o_mu_mu_sigma_nu, Named("mu_mu_nu_nu") = o_mu_mu_nu_nu, Named("mu_sigma_sigma_sigma") = o_mu_sigma_sigma_sigma, Named("mu_sigma_sigma_nu") = o_mu_sigma_sigma_nu, Named("mu_sigma_nu_nu") = o_mu_sigma_nu_nu, Named("mu_nu_nu_nu") = o_mu_nu_nu_nu, Named("sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma, Named("sigma_sigma_sigma_nu") = o_sigma_sigma_sigma_nu, Named("sigma_sigma_nu_nu") = o_sigma_sigma_nu_nu, Named("sigma_nu_nu_nu") = o_sigma_nu_nu_nu, Named("nu_nu_nu_nu") = o_nu_nu_nu_nu);
}

// order 5 of log f in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t2_deriv5_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
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
  struct Prm { double m, s, ik, V7, w1, w6, w8, w11, w24, w25, w33, w35, w36, w38, w45, w46, w48, w60, w61, w62, w64; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 2.0);
    const double V7 = t2_V7(v);
    const double w1 = std::pow(ik, 2);
    const double w6 = 15*ik;
    const double w8 = 30*ik;
    const double w11 = 24/std::pow(s, 5);
    const double w24 = std::pow(s, -4);
    const double w25 = 6*w24;
    const double w33 = std::pow(ik, 3);
    const double w35 = std::pow(s, -3);
    const double w36 = 12*w35;
    const double w38 = 3*ik;
    const double w45 = std::pow(s, -2);
    const double w46 = w33*w45;
    const double w48 = std::pow(ik, 4);
    const double w60 = std::pow(ik, 5);
    const double w61 = 180*w60;
    const double w62 = 30*w48;
    const double w64 = 360*w60;
    return Prm{m, s, ik, V7, w1, w6, w8, w11, w24, w25, w33, w35, w36, w38, w45, w46, w48, w60, w61, w62, w64};
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
    const double w6 = P.w6; (void) w6;
    const double w8 = P.w8; (void) w8;
    const double w11 = P.w11; (void) w11;
    const double w24 = P.w24; (void) w24;
    const double w25 = P.w25; (void) w25;
    const double w33 = P.w33; (void) w33;
    const double w35 = P.w35; (void) w35;
    const double w36 = P.w36; (void) w36;
    const double w38 = P.w38; (void) w38;
    const double w45 = P.w45; (void) w45;
    const double w46 = P.w46; (void) w46;
    const double w48 = P.w48; (void) w48;
    const double w60 = P.w60; (void) w60;
    const double w61 = P.w61; (void) w61;
    const double w62 = P.w62; (void) w62;
    const double w64 = P.w64; (void) w64;
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
    const double w4 = 3*w3;
    const double w5 = ik*w4;
    const double w7 = t*u;
    const double w9 = -w7*w8 - 10*w7;
    const double w10 = w0*w6 + 5*w0 + w3 + w5 + w9;
    const double w12 = w11*zt;
    const double w13 = w10*w12;
    const double w14 = 5*w3;
    const double w15 = ik*w0;
    const double w16 = 3*w15;
    const double w17 = w14 + w16;
    const double w18 = w0 + w17 + w3*w6 + w9;
    const double w19 = std::pow(t, 3);
    const double w20 = ik*w19;
    const double w21 = std::pow(u, 3);
    const double w22 = -w21;
    const double w23 = t*w3;
    const double w26 = -w6*w7 - 5*w7;
    const double w27 = 2*w3;
    const double w28 = w0 + w27;
    const double w29 = 6*w15;
    const double w30 = w29 + w5;
    const double w31 = 24*zt;
    const double w32 = w24*w31;
    const double w34 = 9*w15 + w26;
    const double w37 = w0*zt;
    const double w39 = u*w0;
    const double w40 = 10*w39;
    const double w41 = u*w15;
    const double w42 = 45*w41;
    const double w43 = ik*w23;
    const double w44 = ik*w7;
    const double w47 = w26 + w4;
    const double w49 = (t*w38 - u)/s;
    const double w50 = std::pow(u, 5);
    const double w51 = std::pow(u, 4);
    const double w52 = t*w51;
    const double w53 = std::pow(t, 4);
    const double w54 = u*w53;
    const double w55 = z*zt;
    const double w56 = w0*w21;
    const double w57 = 5*t*w21;
    const double w58 = u*w19;
    const double w59 = w0*w3;
    const double w63 = w19*w3;
    o_mu_mu_mu_mu_mu[i] = w13*w2;
    o_mu_mu_mu_mu_sigma[i] = -w11*w18*w20;
    o_mu_mu_mu_mu_nu[i] = w2*w25*(60*ik*u*w0 + 15*u*w0 - w19 - 6*w20 - w22 - w23*w8 - 15*w23);
    o_mu_mu_mu_sigma_sigma[i] = -w10*w12*w15;
    o_mu_mu_mu_sigma_nu[i] = w2*w32*(-w26 - w28 - w30);
    o_mu_mu_mu_nu_nu[i] = w33*w36*w37*(-w28 - w34);
    o_mu_mu_sigma_sigma_sigma[i] = w11*w18*w19;
    o_mu_mu_sigma_sigma_nu[i] = w24*w29*(45*ik*t*w3 + 3*ik*w19 + 20*t*w3 - w21*w38 - 2*w21 - w40 - w42);
    o_mu_mu_sigma_nu_nu[i] = w2*w36*(t*w14 + w20 + w22 - 2*w39 - 13*w41 + 10*w43);
    o_mu_mu_nu_nu_nu[i] = 6*w19*w46*(w17 - 21*w44 - 3*w7);
    o_mu_sigma_sigma_sigma_sigma[i] = w0*w13;
    o_mu_sigma_sigma_sigma_nu[i] = w15*w32*(6*ik*w3 + w16 + w47);
    o_mu_sigma_sigma_nu_nu[i] = 4*t*w1*w35*zt*(9*w20 + w22 + 13*w23 - w40 - w42 + 18*w43);
    o_mu_sigma_nu_nu_nu[i] = 12*w37*w46*(w16 + w27 - 9*w44 - 2*w7);
    o_mu_nu_nu_nu_nu[i] = w19*w31*w48*w49;
    o_sigma_sigma_sigma_sigma_sigma[i] = w11*(-std::pow(t, 5) + w14*w19 + 2*w50 + w51*w55 + 10*w52 + 15*w53*w55 + 40*w54 + w55*w57 + 5*w55*w58 + 10*w55*w59 + 20*w56);
    o_sigma_sigma_sigma_sigma_nu[i] = u*w25*(6*ik*w0*w3 + 30*ik*w53 - 60*u*w20 - w51 - w57 - 35*w58 - 7*w59);
    o_sigma_sigma_sigma_nu_nu[i] = w36*w41*(w30 + w47);
    o_sigma_sigma_nu_nu_nu[i] = 6*w1*w39*w45*(w34 + w4);
    o_sigma_nu_nu_nu_nu[i] = 24*w33*w49*w58;
    o_nu_nu_nu_nu_nu[i] = V7 - 3*w48*w50 - 15*w48*w52 + 36*w50*w60 + w52*w61 + w54*w61 - w56*w62 + w56*w64 - w62*w63 + w63*w64;
  });
  return List::create(Named("mu_mu_mu_mu_mu") = o_mu_mu_mu_mu_mu, Named("mu_mu_mu_mu_sigma") = o_mu_mu_mu_mu_sigma, Named("mu_mu_mu_mu_nu") = o_mu_mu_mu_mu_nu, Named("mu_mu_mu_sigma_sigma") = o_mu_mu_mu_sigma_sigma, Named("mu_mu_mu_sigma_nu") = o_mu_mu_mu_sigma_nu, Named("mu_mu_mu_nu_nu") = o_mu_mu_mu_nu_nu, Named("mu_mu_sigma_sigma_sigma") = o_mu_mu_sigma_sigma_sigma, Named("mu_mu_sigma_sigma_nu") = o_mu_mu_sigma_sigma_nu, Named("mu_mu_sigma_nu_nu") = o_mu_mu_sigma_nu_nu, Named("mu_mu_nu_nu_nu") = o_mu_mu_nu_nu_nu, Named("mu_sigma_sigma_sigma_sigma") = o_mu_sigma_sigma_sigma_sigma, Named("mu_sigma_sigma_sigma_nu") = o_mu_sigma_sigma_sigma_nu, Named("mu_sigma_sigma_nu_nu") = o_mu_sigma_sigma_nu_nu, Named("mu_sigma_nu_nu_nu") = o_mu_sigma_nu_nu_nu, Named("mu_nu_nu_nu_nu") = o_mu_nu_nu_nu_nu, Named("sigma_sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma_sigma, Named("sigma_sigma_sigma_sigma_nu") = o_sigma_sigma_sigma_sigma_nu, Named("sigma_sigma_sigma_nu_nu") = o_sigma_sigma_sigma_nu_nu, Named("sigma_sigma_nu_nu_nu") = o_sigma_sigma_nu_nu_nu, Named("sigma_nu_nu_nu_nu") = o_sigma_nu_nu_nu_nu, Named("nu_nu_nu_nu_nu") = o_nu_nu_nu_nu_nu);
}

// first derivatives of the expected Hessian
// [[Rcpp::export]]
List student_t2_dexpected1_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
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
  struct Prm { double m, s, ik, V8, w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 2.0);
    const double V8 = t2_V8(v);
    const double w0 = std::pow(ik, 2);
    const double w1 = 6*w0;
    const double w2 = 5*ik;
    const double w3 = w2 + 1;
    const double w4 = 1/(std::pow(s, 3)*w3);
    const double w5 = std::pow(s, -2);
    const double w6 = w5/(10*ik + 25*w0 + 1);
    const double w7 = std::pow(ik, 3);
    const double w8 = 6*w7;
    const double w9 = 15*w0;
    const double w10 = std::pow(ik, 4);
    const double w11 = 16*ik;
    return Prm{m, s, ik, V8, w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11};
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
    const double w5 = P.w5; (void) w5;
    const double w6 = P.w6; (void) w6;
    const double w7 = P.w7; (void) w7;
    const double w8 = P.w8; (void) w8;
    const double w9 = P.w9; (void) w9;
    const double w10 = P.w10; (void) w10;
    const double w11 = P.w11; (void) w11;
    o_mu_mu_mu[i] = 0;
    o_mu_mu_sigma[i] = 2*w4*(w1 + w3);
    o_mu_mu_nu[i] = w6*w8*(w2 + 2);
    o_sigma_sigma_mu[i] = 0;
    o_sigma_sigma_sigma[i] = 4*w4*(2*ik + 1);
    o_sigma_sigma_nu[i] = -w1*w6;
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
    o_sigma_nu_sigma[i] = w5*w8/(8*ik + w9 + 1);
    o_sigma_nu_nu[i] = 6*w10*(w11 + w9 + 3)/(s*(94*w0 + 225*w10 + w11 + 240*w7 + 1));
  });
  return List::create(Named("mu_mu_mu") = o_mu_mu_mu, Named("mu_mu_sigma") = o_mu_mu_sigma, Named("mu_mu_nu") = o_mu_mu_nu, Named("sigma_sigma_mu") = o_sigma_sigma_mu, Named("sigma_sigma_sigma") = o_sigma_sigma_sigma, Named("sigma_sigma_nu") = o_sigma_sigma_nu, Named("nu_nu_mu") = o_nu_nu_mu, Named("nu_nu_sigma") = o_nu_nu_sigma, Named("nu_nu_nu") = o_nu_nu_nu, Named("mu_sigma_mu") = o_mu_sigma_mu, Named("mu_sigma_sigma") = o_mu_sigma_sigma, Named("mu_sigma_nu") = o_mu_sigma_nu, Named("mu_nu_mu") = o_mu_nu_mu, Named("mu_nu_sigma") = o_mu_nu_sigma, Named("mu_nu_nu") = o_mu_nu_nu, Named("sigma_nu_mu") = o_sigma_nu_mu, Named("sigma_nu_sigma") = o_sigma_nu_sigma, Named("sigma_nu_nu") = o_sigma_nu_nu);
}

// second derivatives of the expected Hessian
// [[Rcpp::export]]
List student_t2_dexpected2_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
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
  struct Prm { double m, s, ik, V9, w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14, w15, w16, w17; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 2.0);
    const double V9 = t2_V9(v);
    const double w0 = std::pow(ik, 2);
    const double w1 = 5*ik;
    const double w2 = w1 + 1;
    const double w3 = 1/(std::pow(s, 4)*w2);
    const double w4 = std::pow(ik, 4);
    const double w5 = 15*ik;
    const double w6 = 25*w0;
    const double w7 = std::pow(s, -2);
    const double w8 = std::pow(ik, 3);
    const double w9 = w7/(75*w0 + w5 + 125*w8 + 1);
    const double w10 = 12*w8;
    const double w11 = std::pow(s, -3);
    const double w12 = w11/(10*ik + w6 + 1);
    const double w13 = 15*w0;
    const double w14 = std::pow(ik, 5);
    const double w15 = 237*w0;
    const double w16 = 225*w4;
    const double w17 = 16*ik;
    return Prm{m, s, ik, V9, w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14, w15, w16, w17};
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
    const double w9 = P.w9; (void) w9;
    const double w10 = P.w10; (void) w10;
    const double w11 = P.w11; (void) w11;
    const double w12 = P.w12; (void) w12;
    const double w13 = P.w13; (void) w13;
    const double w14 = P.w14; (void) w14;
    const double w15 = P.w15; (void) w15;
    const double w16 = P.w16; (void) w16;
    const double w17 = P.w17; (void) w17;
    o_mu_mu_mu_mu[i] = 0;
    o_mu_mu_sigma_sigma[i] = -6*w3*(6*w0 + w2);
    o_mu_mu_nu_nu[i] = -12*w4*w9*(w5 + w6 + 3);
    o_mu_mu_mu_sigma[i] = 0;
    o_mu_mu_mu_nu[i] = 0;
    o_mu_mu_sigma_nu[i] = -w10*w12*(w1 + 2);
    o_sigma_sigma_mu_mu[i] = 0;
    o_sigma_sigma_sigma_sigma[i] = -12*w3*(2*ik + 1);
    o_sigma_sigma_nu_nu[i] = w10*w9;
    o_sigma_sigma_mu_sigma[i] = 0;
    o_sigma_sigma_mu_nu[i] = 0;
    o_sigma_sigma_sigma_nu[i] = 12*w0*w12;
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
    o_sigma_nu_sigma_sigma[i] = -w10*w11/(8*ik + w13 + 1);
    o_sigma_nu_nu_nu[i] = -12*w14*(64*ik + w15 + w16 + 360*w8 + 6)/(s*(3375*std::pow(ik, 6) + 24*ik + 5400*w14 + w15 + 3555*w4 + 1232*w8 + 1));
    o_sigma_nu_mu_sigma[i] = 0;
    o_sigma_nu_mu_nu[i] = 0;
    o_sigma_nu_sigma_nu[i] = -6*w4*w7*(w13 + w17 + 3)/(94*w0 + w16 + w17 + 240*w8 + 1);
  });
  return List::create(Named("mu_mu_mu_mu") = o_mu_mu_mu_mu, Named("mu_mu_sigma_sigma") = o_mu_mu_sigma_sigma, Named("mu_mu_nu_nu") = o_mu_mu_nu_nu, Named("mu_mu_mu_sigma") = o_mu_mu_mu_sigma, Named("mu_mu_mu_nu") = o_mu_mu_mu_nu, Named("mu_mu_sigma_nu") = o_mu_mu_sigma_nu, Named("sigma_sigma_mu_mu") = o_sigma_sigma_mu_mu, Named("sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma, Named("sigma_sigma_nu_nu") = o_sigma_sigma_nu_nu, Named("sigma_sigma_mu_sigma") = o_sigma_sigma_mu_sigma, Named("sigma_sigma_mu_nu") = o_sigma_sigma_mu_nu, Named("sigma_sigma_sigma_nu") = o_sigma_sigma_sigma_nu, Named("nu_nu_mu_mu") = o_nu_nu_mu_mu, Named("nu_nu_sigma_sigma") = o_nu_nu_sigma_sigma, Named("nu_nu_nu_nu") = o_nu_nu_nu_nu, Named("nu_nu_mu_sigma") = o_nu_nu_mu_sigma, Named("nu_nu_mu_nu") = o_nu_nu_mu_nu, Named("nu_nu_sigma_nu") = o_nu_nu_sigma_nu, Named("mu_sigma_mu_mu") = o_mu_sigma_mu_mu, Named("mu_sigma_sigma_sigma") = o_mu_sigma_sigma_sigma, Named("mu_sigma_nu_nu") = o_mu_sigma_nu_nu, Named("mu_sigma_mu_sigma") = o_mu_sigma_mu_sigma, Named("mu_sigma_mu_nu") = o_mu_sigma_mu_nu, Named("mu_sigma_sigma_nu") = o_mu_sigma_sigma_nu, Named("mu_nu_mu_mu") = o_mu_nu_mu_mu, Named("mu_nu_sigma_sigma") = o_mu_nu_sigma_sigma, Named("mu_nu_nu_nu") = o_mu_nu_nu_nu, Named("mu_nu_mu_sigma") = o_mu_nu_mu_sigma, Named("mu_nu_mu_nu") = o_mu_nu_mu_nu, Named("mu_nu_sigma_nu") = o_mu_nu_sigma_nu, Named("sigma_nu_mu_mu") = o_sigma_nu_mu_mu, Named("sigma_nu_sigma_sigma") = o_sigma_nu_sigma_sigma, Named("sigma_nu_nu_nu") = o_sigma_nu_nu_nu, Named("sigma_nu_mu_sigma") = o_sigma_nu_mu_sigma, Named("sigma_nu_mu_nu") = o_sigma_nu_mu_nu, Named("sigma_nu_sigma_nu") = o_sigma_nu_sigma_nu);
}

// 1 in y and 1 in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t2_cross_y_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu(n);
  NumericVector o_sigma(n);
  NumericVector o_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, w0; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 2.0);
    const double w0 = 3*ik;
    return Prm{m, s, ik, w0};
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
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    // u = q/(1 + q) in [0, 1], and zt = z t, where q itself may overflow
    const double u = (q < 1.0) ? q * t : 1.0 / (1.0 + 1.0 / q);
    const double zt = (std::fabs(z) <= 1.0) ? z * t : 1.0 / (1.0 / z + z * ik);
    (void) t; (void) q; (void) u; (void) zt;
    const double w1 = t*w0 - u;
    const double w2 = t/std::pow(s, 2);
    o_mu[i] = w2*(t - u*w0 + w1);
    o_sigma[i] = 2*w2*zt*(w0 + 1);
    o_nu[i] = ik*w1*zt/s;
  });
  return List::create(Named("mu") = o_mu, Named("sigma") = o_sigma, Named("nu") = o_nu);
}

// 2 in y and 1 in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t2_cross2_y_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu(n);
  NumericVector o_sigma(n);
  NumericVector o_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, w1, w3; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 2.0);
    const double w1 = 9*ik;
    const double w3 = 2/std::pow(s, 3);
    return Prm{m, s, ik, w1, w3};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double w1 = P.w1; (void) w1;
    const double w3 = P.w3; (void) w3;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    // u = q/(1 + q) in [0, 1], and zt = z t, where q itself may overflow
    const double u = (q < 1.0) ? q * t : 1.0 / (1.0 + 1.0 / q);
    const double zt = (std::fabs(z) <= 1.0) ? z * t : 1.0 / (1.0 / z + z * ik);
    (void) t; (void) q; (void) u; (void) zt;
    const double w0 = 3*t;
    const double w2 = t*w1;
    const double w4 = ik*t;
    const double w5 = std::pow(t, 2);
    o_mu[i] = w3*w4*zt*(3*ik*u + u - w0 - w2);
    o_sigma[i] = w3*w5*(ik*w0 + t - u*w1 - 3*u);
    o_nu[i] = w4*(3*ik*w5 + std::pow(u, 2) - u*w0 - u*w2)/std::pow(s, 2);
  });
  return List::create(Named("mu") = o_mu, Named("sigma") = o_sigma, Named("nu") = o_nu);
}

// 1 in y and 2 in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t2_grad_y_hess_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
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
    const double ik = 1.0 / (v - 2.0);
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
    const double w2 = 3*ik*u;
    const double w3 = 9*ik*t + w0 + w1 - w2;
    const double w4 = ik*t;
    const double w6 = w5*zt;
    const double w7 = ik*w0;
    const double w8 = 2*zt;
    const double w9 = std::pow(t, 2);
    const double w10 = w4/std::pow(s, 2);
    o_mu_mu[i] = w3*w4*w6;
    o_sigma_sigma[i] = -t*w3*w6;
    o_nu_nu[i] = std::pow(ik, 2)*t*w8*(-w1 - w7)/s;
    o_mu_sigma[i] = w5*w9*(9*ik*u - t + 3*u - w7);
    o_mu_nu[i] = w10*(9*ik*t*u - 3*ik*w9 + 3*t*u - std::pow(u, 2));
    o_sigma_nu[i] = w10*w8*(2*u + w2 - w7);
  });
  return List::create(Named("mu_mu") = o_mu_mu, Named("sigma_sigma") = o_sigma_sigma, Named("nu_nu") = o_nu_nu, Named("mu_sigma") = o_mu_sigma, Named("mu_nu") = o_mu_nu, Named("sigma_nu") = o_sigma_nu);
}

// 2 in y and 2 in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t2_hess_y_hess_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu(n);
  NumericVector o_sigma_sigma(n);
  NumericVector o_nu_nu(n);
  NumericVector o_mu_sigma(n);
  NumericVector o_mu_nu(n);
  NumericVector o_sigma_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, w1, w8, w9, w10, w11; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 2.0);
    const double w1 = 3*ik;
    const double w8 = std::pow(s, -4);
    const double w9 = 6*w8;
    const double w10 = std::pow(ik, 2);
    const double w11 = std::pow(s, -3);
    return Prm{m, s, ik, w1, w8, w9, w10, w11};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double w1 = P.w1; (void) w1;
    const double w8 = P.w8; (void) w8;
    const double w9 = P.w9; (void) w9;
    const double w10 = P.w10; (void) w10;
    const double w11 = P.w11; (void) w11;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    // u = q/(1 + q) in [0, 1], and zt = z t, where q itself may overflow
    const double u = (q < 1.0) ? q * t : 1.0 / (1.0 + 1.0 / q);
    const double zt = (std::fabs(z) <= 1.0) ? z * t : 1.0 / (1.0 / z + z * ik);
    (void) t; (void) q; (void) u; (void) zt;
    const double w0 = std::pow(u, 2);
    const double w2 = w0*w1;
    const double w3 = std::pow(t, 2);
    const double w4 = t*u;
    const double w5 = ik*w3;
    const double w6 = -18*ik*w4 + w0;
    const double w7 = w2 + w3 - 6*w4 + 3*w5 + w6;
    o_mu_mu[i] = w5*w7*w9;
    o_sigma_sigma[i] = -w3*w7*w9;
    o_nu_nu[i] = 6*w10*w3*(5*ik*t*u + t*u - w0 - w5)/std::pow(s, 2);
    o_mu_sigma[i] = 24*w5*w8*zt*(t*w1 + t - u*w1 - u);
    o_mu_nu[i] = 2*t*w10*w11*zt*(3*w3 - 8*w4 + 18*w5 + w6);
    o_sigma_nu[i] = 6*w11*w5*(8*ik*t*u + 2*t*u - 2*w0 - w2 - w5);
  });
  return List::create(Named("mu_mu") = o_mu_mu, Named("sigma_sigma") = o_sigma_sigma, Named("nu_nu") = o_nu_nu, Named("mu_sigma") = o_mu_sigma, Named("mu_nu") = o_mu_nu, Named("sigma_nu") = o_sigma_nu);
}

// order 1 of log f in y
// [[Rcpp::export]]
List student_t2_dy1_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_y(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 2.0);
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
    o_y[i] = -zt*(3*ik + 1)/s;
  });
  return List::create(Named("y") = o_y);
}

// order 2 of log f in y
// [[Rcpp::export]]
List student_t2_dy2_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_y(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 2.0);
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
    o_y[i] = t*(-3*ik*t + 3*ik*u - t + u)/std::pow(s, 2);
  });
  return List::create(Named("y") = o_y);
}

// order 3 of log f in y
// [[Rcpp::export]]
List student_t2_dy3_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_y(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 2.0);
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
    const double w0 = ik*t;
    o_y[i] = 2*w0*zt*(-3*ik*u + 3*t - u + 9*w0)/std::pow(s, 3);
  });
  return List::create(Named("y") = o_y);
}

// order 4 of log f in y
// [[Rcpp::export]]
List student_t2_dy4_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_y(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, w3; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 2.0);
    const double w3 = 3*ik;
    return Prm{m, s, ik, w3};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double w3 = P.w3; (void) w3;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    // u = q/(1 + q) in [0, 1], and zt = z t, where q itself may overflow
    const double u = (q < 1.0) ? q * t : 1.0 / (1.0 + 1.0 / q);
    const double zt = (std::fabs(z) <= 1.0) ? z * t : 1.0 / (1.0 / z + z * ik);
    (void) t; (void) q; (void) u; (void) zt;
    const double w0 = std::pow(t, 2);
    const double w1 = std::pow(u, 2);
    const double w2 = t*u;
    o_y[i] = 6*ik*w0*(-18*ik*w2 + w0*w3 + w0 + w1*w3 + w1 - 6*w2)/std::pow(s, 4);
  });
  return List::create(Named("y") = o_y);
}
