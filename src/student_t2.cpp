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
    const double DQ = t2_D(q);
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
  struct Prm { double m, s, ik, V1, w0, w6; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 2.0);
    const double V1 = t2_V1(v);
    const double w0 = 3*ik;
    const double w6 = std::pow(ik, 2);
    return Prm{m, s, ik, V1, w0, w6};
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
    const double w6 = P.w6; (void) w6;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    (void) t; (void) q;
    const double w1 = std::pow(t, 2);
    const double w2 = w1/std::pow(s, 2);
    const double w7 = q*w1;
    const double w8 = (q - w0)/s;
    o_mu_mu[i] = d7::student_t2_hess_mu_mu(z, s, ik);
    o_sigma_sigma[i] = d7::student_t2_hess_sigma_sigma(z, s, ik);
    o_nu_nu[i] = d7::student_t2_hess_nu_nu(z, ik, V1);
    o_mu_sigma[i] = -2*w2*z*(w0 + 1);
    o_mu_nu[i] = ik*w1*w8*z;
    o_sigma_nu[i] = w7*w8;
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
  struct Prm { double m, s, ik, V3, w2, w9, w12, w16; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 2.0);
    const double V3 = t2_V3(v);
    const double w2 = 9*ik;
    const double w9 = -3*ik;
    const double w12 = std::pow(s, -2);
    const double w16 = std::pow(ik, 2);
    return Prm{m, s, ik, V3, w2, w9, w12, w16};
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
    const double w9 = P.w9; (void) w9;
    const double w12 = P.w12; (void) w12;
    const double w16 = P.w16; (void) w16;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    (void) t; (void) q;
    const double w0 = 3*q;
    const double w1 = ik*w0;
    const double w3 = w1 - w2;
    const double w4 = q + w3 - 3;
    const double w5 = std::pow(t, 3);
    const double w6 = 2*w5;
    const double w7 = w6/std::pow(s, 3);
    const double w8 = w7*z;
    const double w10 = q*w2 + w0 + w9;
    const double w11 = std::pow(q, 2);
    const double w13 = w12*w5;
    const double w14 = 2*q;
    const double w15 = w6*z;
    const double w17 = (3*ik - q)/s;
    const double w18 = std::pow(q, 3);
    const double w19 = std::pow(z, 2);
    const double w20 = std::pow(ik, 3)*w5;
    const double w21 = 9*w20;
    const double w22 = w16*w5;
    o_mu_mu_mu[i] = ik*w4*w8;
    o_mu_mu_sigma[i] = w7*(1 - w10);
    o_mu_mu_nu[i] = ik*w13*(-w10 + w11);
    o_mu_sigma_sigma[i] = -w4*w8;
    o_mu_sigma_nu[i] = ik*w12*w15*(-w1 - w14 - w9);
    o_mu_nu_nu[i] = w15*w16*w17;
    o_sigma_sigma_sigma[i] = w7*(15*q + w0*w19 + w11*w19 + 6*w11 + 2*w18 + 6*w19 - 1);
    o_sigma_sigma_nu[i] = q*w13*(-5*q - w11 - w3);
    o_sigma_nu_nu[i] = ik*w14*w17*w5;
    o_nu_nu_nu[i] = V3 + q*w21 + w11*w21 - 3.0/2.0*w11*w22 + 3*w18*w20 - 1.0/2.0*w18*w22;
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
  struct Prm { double m, s, ik, V5, w0, w4, w8, w12, w16, w17, w20, w23, w24, w25, w28, w36; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 2.0);
    const double V5 = t2_V5(v);
    const double w0 = 3*ik;
    const double w4 = 18*ik;
    const double w8 = std::pow(s, -4);
    const double w12 = -w0;
    const double w16 = std::pow(s, -3);
    const double w17 = std::pow(ik, 2);
    const double w20 = ik*w16;
    const double w23 = std::pow(s, -2);
    const double w24 = w17*w23;
    const double w25 = 9*ik;
    const double w28 = std::pow(ik, 3);
    const double w36 = std::pow(ik, 4);
    return Prm{m, s, ik, V5, w0, w4, w8, w12, w16, w17, w20, w23, w24, w25, w28, w36};
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
    const double w4 = P.w4; (void) w4;
    const double w8 = P.w8; (void) w8;
    const double w12 = P.w12; (void) w12;
    const double w16 = P.w16; (void) w16;
    const double w17 = P.w17; (void) w17;
    const double w20 = P.w20; (void) w20;
    const double w23 = P.w23; (void) w23;
    const double w24 = P.w24; (void) w24;
    const double w25 = P.w25; (void) w25;
    const double w28 = P.w28; (void) w28;
    const double w36 = P.w36; (void) w36;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    (void) t; (void) q;
    const double w1 = std::pow(q, 2);
    const double w2 = w0*w1;
    const double w3 = 6*q;
    const double w5 = q*w4;
    const double w6 = w1 - w5;
    const double w7 = w0 + w2 - w3 + w6 + 1;
    const double w9 = std::pow(t, 4);
    const double w10 = 6*w9;
    const double w11 = w10*w8;
    const double w13 = q*w0 + q + w12 - 1;
    const double w14 = w9*z;
    const double w15 = w14*w8;
    const double w18 = 2*w14;
    const double w19 = 2*w1 + w2;
    const double w21 = 5*q;
    const double w22 = -w1;
    const double w26 = 2*q;
    const double w27 = 4*w9;
    const double w29 = w10*w28;
    const double w30 = (q - w0)/s;
    const double w31 = std::pow(q, 3);
    const double w32 = std::pow(q, 4);
    const double w33 = std::pow(z, 2);
    const double w34 = 4*w1;
    const double w35 = w26*w9;
    const double w37 = w36*w9;
    const double w38 = 36*w37;
    const double w39 = w32*w9;
    o_mu_mu_mu_mu[i] = ik*w11*w7;
    o_mu_mu_mu_sigma[i] = -24*ik*w13*w15;
    o_mu_mu_mu_nu[i] = w16*w17*w18*(-8*q + w4 + w6 + 3);
    o_mu_mu_sigma_sigma[i] = -w11*w7;
    o_mu_mu_sigma_nu[i] = w10*w20*(8*ik*q - ik + 2*q - w19);
    o_mu_mu_nu_nu[i] = w10*w24*(ik*w21 - ik + q + w22);
    o_mu_sigma_sigma_sigma[i] = 24*w13*w15;
    o_mu_sigma_sigma_nu[i] = w18*w20*(24*ik*q + 10*q - w19 - w25);
    o_mu_sigma_nu_nu[i] = w24*w27*z*(ik*w3 + w12 + w22 + w26);
    o_mu_nu_nu_nu[i] = w29*w30*z;
    o_sigma_sigma_sigma_sigma[i] = w11*(-26*q - 9*w1 - w21*w33 - w31*w33 - 8*w31 - 2*w32 - w33*w34 - 10*w33 + 1);
    o_sigma_sigma_sigma_nu[i] = w16*w35*(15*q + w31 + w34 - w4 + w5);
    o_sigma_sigma_nu_nu[i] = ik*w23*w35*(q*w25 + w21 + w22 - w25);
    o_sigma_nu_nu_nu[i] = w17*w3*w30*w9;
    o_nu_nu_nu_nu[i] = V5 - q*w38 + w1*w29 - 54*w1*w37 + w27*w28*w31 + w28*w39 - w31*w38 - 9*w36*w39;
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
  struct Prm { double m, s, ik, V7, w0, w2, w4, w7, w17, w22, w33, w34, w35, w38, w42, w43, w46; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 2.0);
    const double V7 = t2_V7(v);
    const double w0 = std::pow(ik, 2);
    const double w2 = 3*ik;
    const double w4 = 15*ik;
    const double w7 = 30*ik;
    const double w17 = 6*ik;
    const double w22 = std::pow(s, -4);
    const double w33 = std::pow(ik, 3);
    const double w34 = 9*ik;
    const double w35 = std::pow(s, -3);
    const double w38 = -w2;
    const double w42 = std::pow(s, -2);
    const double w43 = w33*w42;
    const double w46 = std::pow(ik, 4);
    return Prm{m, s, ik, V7, w0, w2, w4, w7, w17, w22, w33, w34, w35, w38, w42, w43, w46};
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
    const double w0 = P.w0; (void) w0;
    const double w2 = P.w2; (void) w2;
    const double w4 = P.w4; (void) w4;
    const double w7 = P.w7; (void) w7;
    const double w17 = P.w17; (void) w17;
    const double w22 = P.w22; (void) w22;
    const double w33 = P.w33; (void) w33;
    const double w34 = P.w34; (void) w34;
    const double w35 = P.w35; (void) w35;
    const double w38 = P.w38; (void) w38;
    const double w42 = P.w42; (void) w42;
    const double w43 = P.w43; (void) w43;
    const double w46 = P.w46; (void) w46;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    (void) t; (void) q;
    const double w1 = std::pow(q, 2);
    const double w3 = w1*w2;
    const double w5 = 10*q;
    const double w6 = -w5;
    const double w8 = -q*w7 + w6;
    const double w9 = w1 + w3 + w4 + w8 + 5;
    const double w10 = std::pow(t, 5);
    const double w11 = 24*w10;
    const double w12 = w11/std::pow(s, 5);
    const double w13 = w12*w9*z;
    const double w14 = 5*w1;
    const double w15 = w1*w4 + w14 + w2 + w8 + 1;
    const double w16 = ik*w12;
    const double w18 = ik*q;
    const double w19 = 60*w18;
    const double w20 = std::pow(q, 3);
    const double w21 = -w20;
    const double w23 = w10*w22;
    const double w24 = 6*w23;
    const double w25 = 2*w1;
    const double w26 = 5*q;
    const double w27 = q*w4;
    const double w28 = -w26 - w27;
    const double w29 = w25 + w28 + 1;
    const double w30 = w17 + w3;
    const double w31 = w0*z;
    const double w32 = w11*w22;
    const double w36 = 12*w10;
    const double w37 = w35*w36;
    const double w39 = 45*w18;
    const double w40 = 2*q;
    const double w41 = ik*w1;
    const double w44 = 6*w10;
    const double w45 = 3*w1 + w28;
    const double w47 = w11*(-q + w2)/s;
    const double w48 = std::pow(q, 4);
    const double w49 = std::pow(q, 5);
    const double w50 = std::pow(z, 2);
    const double w51 = 5*w20;
    const double w52 = std::pow(ik, 5)*w10;
    const double w53 = 180*w52;
    const double w54 = w10*w46;
    const double w55 = 30*w54;
    const double w56 = 360*w52;
    o_mu_mu_mu_mu_mu[i] = w0*w13;
    o_mu_mu_mu_mu_sigma[i] = -w15*w16;
    o_mu_mu_mu_mu_nu[i] = w0*w24*(15*q - w1*w7 - 15*w1 - w17 + w19 - w21 - 1);
    o_mu_mu_mu_sigma_sigma[i] = -w16*w9*z;
    o_mu_mu_mu_sigma_nu[i] = w31*w32*(-w29 - w30);
    o_mu_mu_mu_nu_nu[i] = w33*w37*z*(-w29 - w34);
    o_mu_mu_sigma_sigma_sigma[i] = w12*w15;
    o_mu_mu_sigma_sigma_nu[i] = w17*w23*(45*ik*w1 + 20*w1 - w2*w20 - 2*w20 - w38 - w39 - w5);
    o_mu_mu_sigma_nu_nu[i] = w0*w37*(ik + w14 - 13*w18 + w21 - w40 + 10*w41);
    o_mu_mu_nu_nu_nu[i] = w43*w44*(-3*q + w14 - 21*w18 - w38);
    o_mu_sigma_sigma_sigma_sigma[i] = w13;
    o_mu_sigma_sigma_sigma_nu[i] = ik*w32*z*(w1*w17 + w2 + w45);
    o_mu_sigma_sigma_nu_nu[i] = 4*w10*w31*w35*(13*w1 + w21 + w34 - w39 + 18*w41 + w6);
    o_mu_sigma_nu_nu_nu[i] = w36*w43*z*(-q*w34 + w25 - w38 - w40);
    o_mu_nu_nu_nu_nu[i] = w46*w47*z;
    o_sigma_sigma_sigma_sigma_sigma[i] = w12*(40*q + 10*w1*w50 + w14 + 20*w20 + w26*w50 + w48*w50 + 10*w48 + 2*w49 + w50*w51 + 15*w50 - 1);
    o_sigma_sigma_sigma_sigma_nu[i] = q*w24*(6*ik*w1 + 30*ik - 35*q - 7*w1 - w19 - w48 - w51);
    o_sigma_sigma_sigma_nu_nu[i] = w18*w37*(w30 + w45);
    o_sigma_sigma_nu_nu_nu[i] = q*w0*w42*w44*(3*w1 - w26 - w27 + w34);
    o_sigma_nu_nu_nu_nu[i] = q*w33*w47;
    o_nu_nu_nu_nu_nu[i] = V7 + q*w53 - w1*w55 + w1*w56 - w20*w55 + w20*w56 + w48*w53 - 15*w48*w54 + 36*w49*w52 - 3*w49*w54;
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
    o_mu_mu_mu[i] = d7::student_t2_dexpected_mu_mu_mu();
    o_mu_mu_sigma[i] = 2*w4*(w1 + w3);
    o_mu_mu_nu[i] = w6*w8*(w2 + 2);
    o_sigma_sigma_mu[i] = 0;
    o_sigma_sigma_sigma[i] = d7::student_t2_dexpected_sigma_sigma_sigma(s, ik);
    o_sigma_sigma_nu[i] = -w1*w6;
    o_nu_nu_mu[i] = 0;
    o_nu_nu_sigma[i] = 0;
    o_nu_nu_nu[i] = d7::student_t2_dexpected_nu_nu_nu(V8);
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
    (void) t; (void) q;
    const double w1 = std::pow(t, 2);
    const double w2 = w1/std::pow(s, 2);
    o_mu[i] = w2*(3*ik - q*w0 - q + 1);
    o_sigma[i] = 2*w2*z*(w0 + 1);
    o_nu[i] = ik*w1*z*(-q + w0)/s;
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
  struct Prm { double m, s, ik, w0; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 2.0);
    const double w0 = 9*ik;
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
    (void) t; (void) q;
    const double w1 = 3*q;
    const double w2 = std::pow(t, 3);
    const double w3 = 2*w2/std::pow(s, 3);
    const double w4 = -3*ik + q*w0 + w1;
    o_mu[i] = ik*w3*z*(ik*w1 + q - w0 - 3);
    o_sigma[i] = w3*(1 - w4);
    o_nu[i] = ik*w2*(std::pow(q, 2) - w4)/std::pow(s, 2);
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
  struct Prm { double m, s, ik, w0, w8, w10; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 2.0);
    const double w0 = 3*ik;
    const double w8 = -w0;
    const double w10 = ik/std::pow(s, 2);
    return Prm{m, s, ik, w0, w8, w10};
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
    const double w8 = P.w8; (void) w8;
    const double w10 = P.w10; (void) w10;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    (void) t; (void) q;
    const double w1 = q*w0;
    const double w2 = -9*ik + q + w1 - 3;
    const double w3 = std::pow(t, 3);
    const double w4 = 2*w3;
    const double w5 = w4/std::pow(s, 3);
    const double w6 = w5*z;
    const double w7 = w4*z;
    const double w9 = 9*ik*q + 3*q + w8;
    o_mu_mu[i] = -ik*w2*w6;
    o_sigma_sigma[i] = w2*w6;
    o_nu_nu[i] = std::pow(ik, 2)*w7*(q - w0)/s;
    o_mu_sigma[i] = w5*(w9 - 1);
    o_mu_nu[i] = w10*w3*(-std::pow(q, 2) + w9);
    o_sigma_nu[i] = w10*w7*(2*q + w1 + w8);
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
  struct Prm { double m, s, ik, w1, w3, w6, w10, w12; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 2.0);
    const double w1 = 3*ik;
    const double w3 = 18*ik;
    const double w6 = std::pow(s, -4);
    const double w10 = std::pow(ik, 2);
    const double w12 = std::pow(s, -3);
    return Prm{m, s, ik, w1, w3, w6, w10, w12};
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
    const double w6 = P.w6; (void) w6;
    const double w10 = P.w10; (void) w10;
    const double w12 = P.w12; (void) w12;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    (void) t; (void) q;
    const double w0 = std::pow(q, 2);
    const double w2 = w0*w1;
    const double w4 = -q*w3 + w0;
    const double w5 = -6*q + w1 + w2 + w4 + 1;
    const double w7 = std::pow(t, 4);
    const double w8 = 6*w7;
    const double w9 = w6*w8;
    const double w11 = w7*z;
    o_mu_mu[i] = ik*w5*w9;
    o_sigma_sigma[i] = -w5*w9;
    o_nu_nu[i] = w10*w8*(5*ik*q - ik + q - w0)/std::pow(s, 2);
    o_mu_sigma[i] = 24*ik*w11*w6*(-q*w1 - q + w1 + 1);
    o_mu_nu[i] = 2*w10*w11*w12*(-8*q + w3 + w4 + 3);
    o_sigma_nu[i] = ik*w12*w8*(8*ik*q - ik + 2*q - 2*w0 - w2);
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
    (void) t; (void) q;
    o_y[i] = -t*z*(3*ik + 1)/s;
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
    (void) t; (void) q;
    o_y[i] = std::pow(t, 2)*(q*w0 + q - w0 - 1)/std::pow(s, 2);
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
    (void) t; (void) q;
    o_y[i] = 2*ik*std::pow(t, 3)*z*(-3*ik*q + 9*ik - q + 3)/std::pow(s, 3);
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
    (void) t; (void) q;
    const double w1 = std::pow(q, 2);
    o_y[i] = 6*ik*std::pow(t, 4)*(-18*ik*q - 6*q + w0*w1 + w0 + w1 + 1)/std::pow(s, 4);
  });
  return List::create(Named("y") = o_y);
}
