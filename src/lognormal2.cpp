#include <Rcpp.h>
#include <cmath>
#include "pt_lognormal2.h"
#include "pt_sqrt.h"
using namespace Rcpp;

// The lognormal in the mean m and the variance v of Y. With C = log(y/m)
// and S = log(1 + v/m^2), log y - mu_l = C + S/2 and the log-density is
//   -log y - log(2 pi)/2 - log(S)/2 - (C + S/2)^2 / (2 S).
// The observed derivatives of orders three to five are in lognormal2_obs.cpp.
// Every other component below is a closed form derived offline with sympy
// (stabilita/gen_lognormal2.py), one exported function per order and per
// surface; y enters through C, the variance only through S. Each is
// checked numerically against independent references in the tests.


// order 1 of log f in (mean, var)
// [[Rcpp::export]]
List lognormal2_gradient_cpp(NumericVector y, NumericVector mean, NumericVector var) {
  const int n = y.size();
  const int n_mean = mean.size(), n_var = var.size();
  NumericVector o_mean(n);
  NumericVector o_var(n);
  double m = 0.0, v = 0.0, S = 0.0;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double t3 = 0.0;
  double t4 = 0.0;
  double t6 = 0.0;
  d7::Lognormal2GradientPar P{};
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; v = var[i % n_var];
    P = d7::lognormal2_gradient_par(m, v);
    S = P.S;
    t0 = P.t0;
    t1 = P.t1;
    t2 = P.t2;
    t3 = P.t3;
    t4 = P.t4;
    t6 = P.t6;
  };
  const bool scalar = n_mean == 1 && n_var == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    o_mean[i] = d7::lognormal2_score_mean(P, C);
    o_var[i] = d7::lognormal2_score_var(P, C);
  }
  (void) S;
  return List::create(Named("mean") = o_mean, Named("var") = o_var);
}

// order 2 of log f in (mean, var)
// [[Rcpp::export]]
List lognormal2_hessian_cpp(NumericVector y, NumericVector mean, NumericVector var) {
  const int n = y.size();
  const int n_mean = mean.size(), n_var = var.size();
  NumericVector o_mean_mean(n);
  NumericVector o_var_var(n);
  NumericVector o_mean_var(n);
  double m = 0.0, v = 0.0, S = 0.0;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double t3 = 0.0;
  double t4 = 0.0;
  double t6 = 0.0;
  double t9 = 0.0;
  double t10 = 0.0;
  double t11 = 0.0;
  double t12 = 0.0;
  double t13 = 0.0;
  double t14 = 0.0;
  double t15 = 0.0;
  double t16 = 0.0;
  d7::Lognormal2HessianPar P{};
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; v = var[i % n_var];
    P = d7::lognormal2_hessian_par(m, v);
    S = P.S;
    t0 = P.t0;
    t1 = P.t1;
    t2 = P.t2;
    t3 = P.t3;
    t4 = P.t4;
    t6 = P.t6;
    t9 = P.t9;
    t10 = P.t10;
    t11 = P.t11;
    t12 = P.t12;
    t13 = P.t13;
    t14 = P.t14;
    t15 = P.t15;
    t16 = P.t16;
  };
  const bool scalar = n_mean == 1 && n_var == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t5 = std::pow(C, 2);
    const double t7 = t5*t6;
    o_mean_mean[i] = d7::lognormal2_hess_mean_mean(P, C);
    o_var_var[i] = d7::lognormal2_hess_var_var(P, C);
    o_mean_var[i] = -C*t15 + m*t12*t2 + m*t12*t3 + (1.0/4.0)*m*t12 + 2*t14*t4*t5*t6 - t15 - t16*t2*t5 - 2*t16*t7;
  }
  (void) S;
  return List::create(Named("mean_mean") = o_mean_mean, Named("var_var") = o_var_var, Named("mean_var") = o_mean_var);
}

// order 2, expected: E_ab on the parameter scale (pt_lognormal2.h)
// [[Rcpp::export]]
List lognormal2_expected_hessian_cpp(NumericVector y, NumericVector mean, NumericVector var) {
  const int n = y.size();
  const int n_mean = mean.size(), n_var = var.size();
  NumericVector o_mean_mean(n);
  NumericVector o_var_var(n);
  NumericVector o_mean_var(n);
  double o[3];
  const bool scalar = n_mean == 1 && n_var == 1;
  if (scalar) d7::lognormal2_expected_all(mean[0], var[0], o);
  for (int i = 0; i < n; i++) {
    if (!scalar) d7::lognormal2_expected_all(mean[i % n_mean], var[i % n_var], o);
    o_mean_mean[i] = o[0];
    o_var_var[i] = o[1];
    o_mean_var[i] = o[2];
  }
  return List::create(Named("mean_mean") = o_mean_mean, Named("var_var") = o_var_var, Named("mean_var") = o_mean_var);
}

// order 3, expected: E[L] = mu_l, E[L^2] = mu_l^2 + S
// [[Rcpp::export]]
List lognormal2_deriv3_expected_cpp(NumericVector y, NumericVector mean, NumericVector var) {
  const int n = y.size();
  const int n_mean = mean.size(), n_var = var.size();
  NumericVector o_mean_mean_mean(n);
  NumericVector o_mean_mean_var(n);
  NumericVector o_mean_var_var(n);
  NumericVector o_var_var_var(n);
  double m = 0.0, v = 0.0, S = 0.0;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double t3 = 0.0;
  double t4 = 0.0;
  double t5 = 0.0;
  double t6 = 0.0;
  double t7 = 0.0;
  double t8 = 0.0;
  double t9 = 0.0;
  double t10 = 0.0;
  double t11 = 0.0;
  double t12 = 0.0;
  double t13 = 0.0;
  double t14 = 0.0;
  double t15 = 0.0;
  double t16 = 0.0;
  double t17 = 0.0;
  double t18 = 0.0;
  double t19 = 0.0;
  double t20 = 0.0;
  double t21 = 0.0;
  double t22 = 0.0;
  double t23 = 0.0;
  double t24 = 0.0;
  double t25 = 0.0;
  double t26 = 0.0;
  double t27 = 0.0;
  double t28 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; v = var[i % n_var];
    S = std::log1p(v / (m * m));
    t0 = std::pow(m, 4);
    t1 = std::pow(v, 2);
    t2 = std::pow(m, 2);
    t3 = 2*v;
    t4 = 1.0/(t0 + t1 + t2*t3);
    t5 = 1.0/m;
    t6 = t4*t5;
    t7 = std::pow(m, 3);
    t8 = 1.0/t7;
    t9 = t1*t4*t8;
    t10 = 1.0/S;
    t11 = t10*t4;
    t12 = t11*t5;
    t13 = std::pow(S, -2);
    t14 = std::pow(v, 3);
    t15 = 3*v;
    t16 = 3*t1;
    t17 = 1.0/(std::pow(m, 6) + t0*t15 + t14 + t16*t2);
    t18 = t17*t7;
    t19 = m*t17;
    t20 = t14*t17*t8;
    t21 = t13*t17;
    t22 = t21*v;
    t23 = t17*t5;
    t24 = std::pow(S, -3);
    t25 = 6*t20;
    t26 = t10*t23;
    t27 = 1.0/t2;
    t28 = t10*t17;
  };
  const bool scalar = n_mean == 1 && n_var == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    o_mean_mean_mean[i] = -6*m*t22 - m*t4 + 6*t1*t21*t5 + 18*t1*t26 + 3*t10*t18 + 9*t10*t19*v + t10*t25 + 6*t10*t9 + 12*t12*v - t13*t25 - 12*t13*t9 + t15*t19 + t16*t23 + t18 - 16*t20*t24 + t20 - t3*t6 - t9;
    o_mean_mean_var[i] = 3*t1*t13*t17*t27 + 8*t1*t17*t24*t27 - 1.0/2.0*t1*t27*t28 - 1.0/2.0*t11*t27*v - 5.0/2.0*t11 + t13*t17*t2 + 4*t13*t27*t4*v - 5*t22 - 7.0/2.0*t28*v;
    o_mean_var_var[i] = m*t10*t17 + 2*m*t13*t17 - 1.0/2.0*t12 - t13*t6 - 5.0/2.0*t22*t5 - 4*t23*t24*v - 1.0/2.0*t26*v;
    o_var_var_var[i] = t17*((3.0/4.0)*t10 + (9.0/4.0)*t13 + 2*t24);
  }
  (void) S;
  return List::create(Named("mean_mean_mean") = o_mean_mean_mean, Named("mean_mean_var") = o_mean_mean_var, Named("mean_var_var") = o_mean_var_var, Named("var_var_var") = o_var_var_var);
}

// order 4, expected: E[L] = mu_l, E[L^2] = mu_l^2 + S
// [[Rcpp::export]]
List lognormal2_deriv4_expected_cpp(NumericVector y, NumericVector mean, NumericVector var) {
  const int n = y.size();
  const int n_mean = mean.size(), n_var = var.size();
  NumericVector o_mean_mean_mean_mean(n);
  NumericVector o_mean_mean_mean_var(n);
  NumericVector o_mean_mean_var_var(n);
  NumericVector o_mean_var_var_var(n);
  NumericVector o_var_var_var_var(n);
  double m = 0.0, v = 0.0, S = 0.0;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double t3 = 0.0;
  double t4 = 0.0;
  double t5 = 0.0;
  double t6 = 0.0;
  double t7 = 0.0;
  double t8 = 0.0;
  double t9 = 0.0;
  double t10 = 0.0;
  double t11 = 0.0;
  double t12 = 0.0;
  double t13 = 0.0;
  double t14 = 0.0;
  double t15 = 0.0;
  double t16 = 0.0;
  double t17 = 0.0;
  double t18 = 0.0;
  double t19 = 0.0;
  double t20 = 0.0;
  double t21 = 0.0;
  double t22 = 0.0;
  double t23 = 0.0;
  double t24 = 0.0;
  double t25 = 0.0;
  double t26 = 0.0;
  double t27 = 0.0;
  double t28 = 0.0;
  double t29 = 0.0;
  double t30 = 0.0;
  double t31 = 0.0;
  double t32 = 0.0;
  double t33 = 0.0;
  double t34 = 0.0;
  double t35 = 0.0;
  double t36 = 0.0;
  double t37 = 0.0;
  double t38 = 0.0;
  double t39 = 0.0;
  double t40 = 0.0;
  double t41 = 0.0;
  double t42 = 0.0;
  double t43 = 0.0;
  double t44 = 0.0;
  double t45 = 0.0;
  double t46 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; v = var[i % n_var];
    S = std::log1p(v / (m * m));
    t0 = std::pow(m, 6);
    t1 = std::pow(v, 3);
    t2 = std::pow(m, 4);
    t3 = 3*t2;
    t4 = std::pow(v, 2);
    t5 = std::pow(m, 2);
    t6 = 1.0/(t0 + t1 + t3*v + 3*t4*t5);
    t7 = 1.0/S;
    t8 = t6*t7;
    t9 = 1.0/t2;
    t10 = 1.0/t5;
    t11 = std::pow(S, -3);
    t12 = std::pow(S, -2);
    t13 = 22*t9;
    t14 = std::pow(v, 4);
    t15 = 1.0/(std::pow(m, 8) + 4*t0*v + 4*t1*t5 + t14 + 6*t2*t4);
    t16 = t15*t4;
    t17 = 18*t16;
    t18 = 12*t15;
    t19 = t5*v;
    t20 = t11*t16;
    t21 = t15*t7;
    t22 = t16*t7;
    t23 = t14*t15;
    t24 = t1*t10;
    t25 = std::pow(S, -4);
    t26 = 1.0/m;
    t27 = t26*t6;
    t28 = t26*t8;
    t29 = std::pow(m, 3);
    t30 = 1.0/t29;
    t31 = t30*t4;
    t32 = t11*t6;
    t33 = t12*t6;
    t34 = t11*t18;
    t35 = m*t34;
    t36 = t12*t15;
    t37 = t36*v;
    t38 = t21*v;
    t39 = t1*t30;
    t40 = t15*t39;
    t41 = 10*t12;
    t42 = t12*t16;
    t43 = t10*v;
    t44 = t11*t15;
    t45 = 3*t27;
    t46 = 18*t26;
  };
  const bool scalar = n_mean == 1 && n_var == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    o_mean_mean_mean_mean[i] = 192*t1*t10*t11*t15 + 156*t1*t10*t12*t15 - 96*t1*t11*t6*t9 + 72*t1*t12*t6*t9 - t1*t13*t8 + 3*t1*t6*t9 + 168*t10*t12*t4*t6 + 9*t10*t4*t6 - 64*t10*t4*t8 + 50*t12*t14*t15*t9 + 6*t12*t15*t4 + 60*t12*t15*t5*v - t13*t14*t21 - t15*t3 - t17 - t18*t19 - t18*t24 - 44*t19*t21 - 11*t2*t21 - 48*t20 - 86*t21*t24 - 141*t22 - 144*t23*t25*t9 - 3*t23*t9 + 3*t5*t6 + 9*t6*v - 74*t8*v;
    o_mean_mean_mean_var[i] = 27*m*t37 + 24*m*t38 + 13*m*t8 - 48*t12*t27*v - 96*t20*t26 + t21*t39 + 3*t22*t26 + 72*t25*t40 - 54*t26*t42 + 2*t28*v - 9*t29*t36 + 36*t31*t32 - 12*t31*t33 + t31*t8 + t34*t39 + t35*v - t40*t41;
    o_mean_mean_var_var[i] = -t10*t11*t17 + (1.0/2.0)*t10*t15*t4*t7 - 36*t10*t16*t25 - 3.0/2.0*t10*t42 + (1.0/2.0)*t10*t6*t7*v + 40*t11*t15*v + (55.0/2.0)*t12*t15*v + 9*t12*t6 + (13.0/2.0)*t15*t7*v - 5*t21*t5 - 12*t32*t43 - 3*t33*t43 - 11*t36*t5 - 2*t44*t5 + (9.0/2.0)*t6*t7;
    o_mean_var_var_var[i] = -9.0/2.0*m*t21 - 27.0/2.0*m*t36 + t11*t45 + t12*t45 + t15*t25*t46*v + (13.0/2.0)*t26*t37 + t26*t38 + t28 - t35 + t44*t46*v;
    o_var_var_var_var[i] = t15*(-15*t11 - 9*t25 - t41 - 11.0/4.0*t7);
  }
  (void) S;
  return List::create(Named("mean_mean_mean_mean") = o_mean_mean_mean_mean, Named("mean_mean_mean_var") = o_mean_mean_mean_var, Named("mean_mean_var_var") = o_mean_mean_var_var, Named("mean_var_var_var") = o_mean_var_var_var, Named("var_var_var_var") = o_var_var_var_var);
}

// first derivatives of the expected Hessian (pt_lognormal2.h)
// [[Rcpp::export]]
List lognormal2_dexpected1_cpp(NumericVector y, NumericVector mean, NumericVector var) {
  const int n = y.size();
  const int n_mean = mean.size(), n_var = var.size();
  NumericVector o_mean_mean_mean(n);
  NumericVector o_mean_mean_var(n);
  NumericVector o_var_var_mean(n);
  NumericVector o_var_var_var(n);
  NumericVector o_mean_var_mean(n);
  NumericVector o_mean_var_var(n);
  double o[6];
  const bool scalar = n_mean == 1 && n_var == 1;
  if (scalar) d7::lognormal2_dexpected1_all(mean[0], var[0], o);
  for (int i = 0; i < n; i++) {
    if (!scalar) d7::lognormal2_dexpected1_all(mean[i % n_mean], var[i % n_var], o);
    o_mean_mean_mean[i] = o[0];
    o_mean_mean_var[i] = o[1];
    o_var_var_mean[i] = o[2];
    o_var_var_var[i] = o[3];
    o_mean_var_mean[i] = o[4];
    o_mean_var_var[i] = o[5];
  }
  return List::create(Named("mean_mean_mean") = o_mean_mean_mean, Named("mean_mean_var") = o_mean_mean_var, Named("var_var_mean") = o_var_var_mean, Named("var_var_var") = o_var_var_var, Named("mean_var_mean") = o_mean_var_mean, Named("mean_var_var") = o_mean_var_var);
}

// second derivatives of the expected Hessian (pt_lognormal2.h)
// [[Rcpp::export]]
List lognormal2_dexpected2_cpp(NumericVector y, NumericVector mean, NumericVector var) {
  const int n = y.size();
  const int n_mean = mean.size(), n_var = var.size();
  NumericVector o_mean_mean_mean_mean(n);
  NumericVector o_mean_mean_var_var(n);
  NumericVector o_mean_mean_mean_var(n);
  NumericVector o_var_var_mean_mean(n);
  NumericVector o_var_var_var_var(n);
  NumericVector o_var_var_mean_var(n);
  NumericVector o_mean_var_mean_mean(n);
  NumericVector o_mean_var_var_var(n);
  NumericVector o_mean_var_mean_var(n);
  double o[9];
  const bool scalar = n_mean == 1 && n_var == 1;
  if (scalar) d7::lognormal2_dexpected2_all(mean[0], var[0], o);
  for (int i = 0; i < n; i++) {
    if (!scalar) d7::lognormal2_dexpected2_all(mean[i % n_mean], var[i % n_var], o);
    o_mean_mean_mean_mean[i] = o[0];
    o_mean_mean_var_var[i] = o[1];
    o_mean_mean_mean_var[i] = o[2];
    o_var_var_mean_mean[i] = o[3];
    o_var_var_var_var[i] = o[4];
    o_var_var_mean_var[i] = o[5];
    o_mean_var_mean_mean[i] = o[6];
    o_mean_var_var_var[i] = o[7];
    o_mean_var_mean_var[i] = o[8];
  }
  return List::create(Named("mean_mean_mean_mean") = o_mean_mean_mean_mean, Named("mean_mean_var_var") = o_mean_mean_var_var, Named("mean_mean_mean_var") = o_mean_mean_mean_var, Named("var_var_mean_mean") = o_var_var_mean_mean, Named("var_var_var_var") = o_var_var_var_var, Named("var_var_mean_var") = o_var_var_mean_var, Named("mean_var_mean_mean") = o_mean_var_mean_mean, Named("mean_var_var_var") = o_mean_var_var_var, Named("mean_var_mean_var") = o_mean_var_mean_var);
}

// order 1 of the distribution function in (mean, var)
// [[Rcpp::export]]
List lognormal2_dcdf1_cpp(NumericVector y, NumericVector mean, NumericVector var) {
  const int n = y.size();
  const int n_mean = mean.size(), n_var = var.size();
  NumericVector o_mean(n);
  NumericVector o_var(n);
  double m = 0.0, v = 0.0, S = 0.0;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; v = var[i % n_var];
    S = std::log1p(v / (m * m));
    t0 = 1.0/m;
    t1 = 1.0/S;
    t2 = 1.0/(std::pow(m, 2) + v);
  };
  const bool scalar = n_mean == 1 && n_var == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t3 = 2*C*t1;
    const double t4 = M_SQRT2*std::exp(-1.0/2.0*C)*std::exp(-1.0/8.0*S)*std::exp(-1.0/2.0*std::pow(C, 2)*t1)/(d7::sqrt_cr(M_PI)*d7::sqrt_cr(S));
    o_mean[i] = (1.0/4.0)*t4*(2*C*t0*t1 - m*t2*t3 + m*t2 - 3*t0);
    o_var[i] = (1.0/8.0)*t2*t4*(1 - t3);
  }
  (void) S;
  return List::create(Named("mean") = o_mean, Named("var") = o_var);
}

// order 2 of the distribution function in (mean, var)
// [[Rcpp::export]]
List lognormal2_dcdf2_cpp(NumericVector y, NumericVector mean, NumericVector var) {
  const int n = y.size();
  const int n_mean = mean.size(), n_var = var.size();
  NumericVector o_mean_mean(n);
  NumericVector o_var_var(n);
  NumericVector o_mean_var(n);
  double m = 0.0, v = 0.0, S = 0.0;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double t3 = 0.0;
  double t4 = 0.0;
  double t5 = 0.0;
  double t8 = 0.0;
  double t14 = 0.0;
  double t15 = 0.0;
  double t17 = 0.0;
  double t19 = 0.0;
  double t20 = 0.0;
  double t21 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; v = var[i % n_var];
    S = std::log1p(v / (m * m));
    t0 = std::pow(m, 2);
    t1 = 1.0/t0;
    t2 = 1.0/S;
    t3 = t1*t2;
    t4 = std::pow(S, -2);
    t5 = 1.0/(t0 + v);
    t8 = std::pow(S, -3);
    t14 = 1.0/(std::pow(m, 4) + 2*t0*v + std::pow(v, 2));
    t15 = t0*t14;
    t17 = 4*t2;
    t19 = 1.0/m;
    t20 = t19*t5;
    t21 = m*t14;
  };
  const bool scalar = n_mean == 1 && n_var == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t6 = std::pow(C, 2);
    const double t7 = std::pow(C, 3);
    const double t9 = t7*t8;
    const double t10 = (1.0/2.0)*t9;
    const double t11 = C*t4;
    const double t12 = C*t2*t5;
    const double t13 = t4*t6;
    const double t16 = M_SQRT2*std::exp(-1.0/2.0*C)*std::exp(-1.0/8.0*S)*std::exp(-1.0/2.0*t2*t6)/(d7::sqrt_cr(M_PI)*d7::sqrt_cr(S));
    const double t18 = 8*t9;
    o_mean_mean[i] = t16*((9.0/8.0)*C*t0*t14*t2 + (3.0/2.0)*C*t0*t14*t4 + (3.0/2.0)*C*t1*t4 - 7.0/8.0*C*t3 + (1.0/4.0)*t0*t14*t4*t6 - t1*t10 + (5.0/4.0)*t1*t4*t6 + (3.0/16.0)*t1 - t10*t15 - 3*t11*t5 - 3.0/4.0*t12 - 3.0/2.0*t13*t5 - 1.0/4.0*t15*t2 - 9.0/16.0*t15 + (3.0/2.0)*t2*t5 - 5.0/4.0*t3 + t5*t7*t8 + (5.0/8.0)*t5);
    o_var_var[i] = (1.0/64.0)*t14*t16*(18*C*t2 + 24*C*t4 - t17 - t18 + 4*t4*t6 - 9);
    o_mean_var[i] = (1.0/32.0)*t16*(18*C*m*t14*t2 + 24*C*m*t14*t4 + 4*m*t14*t4*t6 - 24*t11*t20 - 2*t12*t19 - 12*t13*t20 - t17*t21 - t18*t21 + 12*t19*t2*t5 + 8*t19*t5*t7*t8 + 3*t19*t5 - 9*t21);
  }
  (void) S;
  return List::create(Named("mean_mean") = o_mean_mean, Named("var_var") = o_var_var, Named("mean_var") = o_mean_var);
}

// order 3 of the distribution function in (mean, var)
// [[Rcpp::export]]
List lognormal2_dcdf3_cpp(NumericVector y, NumericVector mean, NumericVector var) {
  const int n = y.size();
  const int n_mean = mean.size(), n_var = var.size();
  NumericVector o_mean_mean_mean(n);
  NumericVector o_mean_mean_var(n);
  NumericVector o_mean_var_var(n);
  NumericVector o_var_var_var(n);
  double m = 0.0, v = 0.0, S = 0.0;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double t3 = 0.0;
  double t4 = 0.0;
  double t5 = 0.0;
  double t6 = 0.0;
  double t12 = 0.0;
  double t20 = 0.0;
  double t23 = 0.0;
  double t24 = 0.0;
  double t25 = 0.0;
  double t26 = 0.0;
  double t27 = 0.0;
  double t29 = 0.0;
  double t30 = 0.0;
  double t31 = 0.0;
  double t32 = 0.0;
  double t36 = 0.0;
  double t37 = 0.0;
  double t38 = 0.0;
  double t39 = 0.0;
  double t41 = 0.0;
  double t42 = 0.0;
  double t45 = 0.0;
  double t46 = 0.0;
  double t47 = 0.0;
  double t48 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; v = var[i % n_var];
    S = std::log1p(v / (m * m));
    t0 = std::pow(m, 3);
    t1 = 1.0/t0;
    t2 = std::pow(S, -2);
    t3 = t1*t2;
    t4 = 1.0/S;
    t5 = (21.0/8.0)*t4;
    t6 = std::pow(S, -3);
    t12 = std::pow(S, -4);
    t20 = std::pow(S, -5);
    t23 = 1.0/m;
    t24 = std::pow(m, 2);
    t25 = 1.0/(t24 + v);
    t26 = t23*t25;
    t27 = (45.0/2.0)*t26;
    t29 = std::pow(m, 4);
    t30 = std::pow(v, 2);
    t31 = 1.0/(2*t24*v + t29 + t30);
    t32 = m*t31;
    t36 = (51.0/8.0)*t32;
    t37 = 1.0/(std::pow(m, 6) + 3*t24*t30 + 3*t29*v + std::pow(v, 3));
    t38 = t0*t37;
    t39 = (13.0/8.0)*t38;
    t41 = 1.0/t24;
    t42 = t25*t41;
    t45 = t2*t31;
    t46 = t24*t37;
    t47 = t23*t31;
    t48 = m*t37;
  };
  const bool scalar = n_mean == 1 && n_var == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t7 = C*t6;
    const double t8 = (15.0/2.0)*t7;
    const double t9 = C*t4;
    const double t10 = std::pow(C, 2);
    const double t11 = t10*t6;
    const double t13 = std::pow(C, 3);
    const double t14 = t12*t13;
    const double t15 = 5*t14;
    const double t16 = t13*t6;
    const double t17 = std::pow(C, 4);
    const double t18 = t12*t17;
    const double t19 = std::pow(C, 5);
    const double t21 = t19*t20;
    const double t22 = (1.0/2.0)*t21;
    const double t28 = C*t2;
    const double t33 = t10*t2;
    const double t34 = 15*t14;
    const double t35 = (3.0/2.0)*t21;
    const double t40 = M_SQRT2*std::exp(-1.0/2.0*C)*std::exp(-1.0/8.0*S)*std::exp(-1.0/2.0*t10*t4)/(d7::sqrt_cr(M_PI)*d7::sqrt_cr(S));
    const double t43 = (15.0/4.0)*t7;
    const double t44 = (1.0/4.0)*t21;
    const double t49 = 480*t7;
    const double t50 = 624*t28;
    const double t51 = 306*t9;
    const double t52 = 96*t11;
    const double t53 = 104*t33;
    const double t54 = 32*t21;
    o_mean_mean_mean[i] = t40*(-33.0/4.0*C*t3 + (21.0/2.0)*t1*t11 - t1*t15 + (11.0/4.0)*t1*t16 - 7.0/4.0*t1*t18 + t1*t22 + t1*t5 + t1*t8 + (41.0/32.0)*t1*t9 - 15.0/64.0*t1 - 21.0/8.0*t10*t3 - t11*t27 + (27.0/2.0)*t11*t32 - 3.0/2.0*t11*t38 + t15*t38 - 3.0/4.0*t16*t26 - 21.0/4.0*t16*t32 + (13.0/4.0)*t16*t38 + (15.0/4.0)*t18*t26 - 9.0/4.0*t18*t32 + (1.0/4.0)*t18*t38 + (45.0/4.0)*t2*t26 - 27.0/4.0*t2*t32 + (3.0/4.0)*t2*t38 - t22*t38 + (9.0/4.0)*t26*t28 - 21.0/8.0*t26*t33 + t26*t34 - t26*t35 + t26*t5 + (9.0/32.0)*t26*t9 + (27.0/64.0)*t26 - t27*t7 + (63.0/4.0)*t28*t32 - 39.0/4.0*t28*t38 - 21.0/4.0*t3 - t32*t34 + t32*t35 + (45.0/2.0)*t32*t7 + (135.0/32.0)*t32*t9 - 189.0/64.0*t32 + t33*t36 - t33*t39 - t36*t4 - t38*t8 - 153.0/32.0*t38*t9 + (153.0/64.0)*t38 + t39*t4);
    o_mean_mean_var[i] = t40*((9.0/8.0)*C*t2*t25*t41 + (9.0/2.0)*C*t2*t31 + (7.0/64.0)*C*t25*t4*t41 + (27.0/32.0)*C*t31*t4 + (15.0/2.0)*C*t31*t6 + 2*t10*t2*t31 + (9.0/2.0)*t10*t31*t6 - 15.0/4.0*t11*t42 - 3.0/4.0*t11*t46 + (5.0/2.0)*t12*t13*t24*t37 + (5.0/2.0)*t12*t13*t25*t41 + (1.0/8.0)*t12*t17*t24*t37 + (5.0/8.0)*t12*t17*t25*t41 + (13.0/8.0)*t13*t24*t37*t6 - t15*t31 - 3.0/2.0*t16*t31 - 3.0/8.0*t16*t42 - 3.0/4.0*t18*t31 + (1.0/2.0)*t19*t20*t31 + (3.0/8.0)*t2*t24*t37 + (15.0/8.0)*t2*t25*t41 + (13.0/16.0)*t24*t37*t4 + (153.0/128.0)*t24*t37 + (1.0/16.0)*t25*t4*t41 - 39.0/8.0*t28*t46 - 2*t31*t4 - 45.0/64.0*t31 - 1.0/16.0*t33*t42 - 13.0/16.0*t33*t46 - t42*t43 - t42*t44 - 3.0/128.0*t42 - t43*t46 - t44*t46 - 9.0/4.0*t45 - 153.0/64.0*t46*t9);
    o_mean_var_var[i] = (1.0/256.0)*t40*(240*C*t2*t23*t31 + 18*C*t23*t31*t4 + 480*C*t23*t31*t6 + 320*m*t12*t13*t37 + 16*m*t12*t17*t37 + 208*m*t13*t37*t6 + 48*m*t2*t37 + 104*m*t37*t4 + 153*m*t37 + 120*t10*t2*t23*t31 + 288*t10*t23*t31*t6 - 320*t14*t47 - 80*t16*t47 - 48*t18*t47 + 32*t19*t20*t23*t31 - 144*t23*t45 - 120*t4*t47 - 27*t47 - t48*t49 - t48*t50 - t48*t51 - t48*t52 - t48*t53 - t48*t54);
    o_var_var_var[i] = (1.0/512.0)*t37*t40*(320*t12*t13 + 16*t12*t17 + 208*t13*t6 + 48*t2 + 104*t4 - t49 - t50 - t51 - t52 - t53 - t54 + 153);
  }
  (void) S;
  return List::create(Named("mean_mean_mean") = o_mean_mean_mean, Named("mean_mean_var") = o_mean_mean_var, Named("mean_var_var") = o_mean_var_var, Named("var_var_var") = o_var_var_var);
}

// order 4 of the distribution function in (mean, var)
// [[Rcpp::export]]
List lognormal2_dcdf4_cpp(NumericVector y, NumericVector mean, NumericVector var) {
  const int n = y.size();
  const int n_mean = mean.size(), n_var = var.size();
  NumericVector o_mean_mean_mean_mean(n);
  NumericVector o_mean_mean_mean_var(n);
  NumericVector o_mean_mean_var_var(n);
  NumericVector o_mean_var_var_var(n);
  NumericVector o_var_var_var_var(n);
  double m = 0.0, v = 0.0, S = 0.0;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double t3 = 0.0;
  double t4 = 0.0;
  double t5 = 0.0;
  double t6 = 0.0;
  double t7 = 0.0;
  double t10 = 0.0;
  double t15 = 0.0;
  double t19 = 0.0;
  double t22 = 0.0;
  double t23 = 0.0;
  double t24 = 0.0;
  double t25 = 0.0;
  double t28 = 0.0;
  double t29 = 0.0;
  double t37 = 0.0;
  double t38 = 0.0;
  double t39 = 0.0;
  double t40 = 0.0;
  double t41 = 0.0;
  double t42 = 0.0;
  double t43 = 0.0;
  double t44 = 0.0;
  double t51 = 0.0;
  double t52 = 0.0;
  double t53 = 0.0;
  double t55 = 0.0;
  double t56 = 0.0;
  double t57 = 0.0;
  double t60 = 0.0;
  double t61 = 0.0;
  double t62 = 0.0;
  double t65 = 0.0;
  double t66 = 0.0;
  double t67 = 0.0;
  double t68 = 0.0;
  double t69 = 0.0;
  double t70 = 0.0;
  double t77 = 0.0;
  double t78 = 0.0;
  double t79 = 0.0;
  double t81 = 0.0;
  double t84 = 0.0;
  double t85 = 0.0;
  double t86 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; v = var[i % n_var];
    S = std::log1p(v / (m * m));
    t0 = std::pow(m, 4);
    t1 = 1.0/t0;
    t2 = std::pow(S, -3);
    t3 = t1*t2;
    t4 = std::pow(S, -2);
    t5 = 1.0/S;
    t6 = t1*t5;
    t7 = std::pow(S, -4);
    t10 = std::pow(S, -5);
    t15 = std::pow(S, -6);
    t19 = std::pow(S, -7);
    t22 = std::pow(m, 2);
    t23 = 1.0/t22;
    t24 = 1.0/(t22 + v);
    t25 = t23*t24;
    t28 = std::pow(v, 2);
    t29 = 1.0/(t0 + 2*t22*v + t28);
    t37 = (225.0/2.0)*t29;
    t38 = t29*t4;
    t39 = t29*t5;
    t40 = std::pow(m, 6);
    t41 = std::pow(v, 3);
    t42 = 3*t22;
    t43 = 1.0/(3*t0*v + t28*t42 + t40 + t41);
    t44 = t22*t43;
    t51 = 1.0/(std::pow(m, 8) + 6*t0*t28 + 4*t22*t41 + 4*t40*v + std::pow(v, 4));
    t52 = t0*t51;
    t53 = (15.0/4.0)*t52;
    t55 = std::pow(m, 3);
    t56 = 1.0/t55;
    t57 = t24*t56;
    t60 = 1.0/m;
    t61 = t29*t60;
    t62 = 14400*t61;
    t65 = m*t43;
    t66 = 3825*t51;
    t67 = 960*t2;
    t68 = t51*t55;
    t69 = 2448*t4;
    t70 = 3212*t5;
    t77 = t23*t29;
    t78 = 4800*t77;
    t79 = t2*t43;
    t81 = t4*t43;
    t84 = t22*t51;
    t85 = t43*t60;
    t86 = m*t51;
  };
  const bool scalar = n_mean == 1 && n_var == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t8 = std::pow(C, 2);
    const double t9 = std::pow(C, 3);
    const double t11 = t10*t9;
    const double t12 = std::pow(C, 4);
    const double t13 = t10*t12;
    const double t14 = std::pow(C, 5);
    const double t16 = t10*t14;
    const double t17 = std::pow(C, 6);
    const double t18 = std::pow(C, 7);
    const double t20 = t18*t19;
    const double t21 = (1.0/2.0)*t20;
    const double t26 = C*t7;
    const double t27 = 210*t26;
    const double t30 = t7*t8;
    const double t31 = t7*t9;
    const double t32 = t2*t9;
    const double t33 = t12*t7;
    const double t34 = t14*t15;
    const double t35 = 42*t34;
    const double t36 = t15*t17;
    const double t45 = C*t2;
    const double t46 = C*t4;
    const double t47 = C*t5;
    const double t48 = t2*t8;
    const double t49 = t4*t8;
    const double t50 = t36*t43;
    const double t54 = M_SQRT2*std::exp(-1.0/2.0*C)*std::exp(-1.0/8.0*S)*std::exp(-1.0/2.0*t5*t8)/(d7::sqrt_cr(M_PI)*d7::sqrt_cr(S));
    const double t58 = 13440*t26;
    const double t59 = 2688*t34;
    const double t63 = C*t38;
    const double t64 = C*t39;
    const double t71 = 13440*t11;
    const double t72 = 16320*t31;
    const double t73 = 6424*t32;
    const double t74 = 960*t13;
    const double t75 = 816*t33;
    const double t76 = 128*t20;
    const double t80 = C*t79;
    const double t82 = C*t81;
    const double t83 = t79*t8;
    o_mean_mean_mean_mean[i] = t54*((765.0/8.0)*C*t0*t2*t51 + (2409.0/32.0)*C*t0*t4*t51 + (3825.0/128.0)*C*t0*t5*t51 + (105.0/2.0)*C*t0*t51*t7 + (969.0/32.0)*C*t1*t4 + (105.0/2.0)*C*t1*t7 + (225.0/2.0)*C*t2*t23*t24 + (495.0/4.0)*C*t2*t29 + (51.0/8.0)*C*t23*t24*t4 + (1.0/32.0)*C*t23*t24*t5 + (243.0/16.0)*C*t29*t4 + (243.0/64.0)*C*t29*t5 + 315*C*t29*t7 - 675.0/8.0*C*t3 - 399.0/128.0*C*t6 + (51.0/8.0)*t0*t10*t14*t51 + (21.0/2.0)*t0*t14*t15*t51 + (1.0/4.0)*t0*t15*t17*t51 + (153.0/8.0)*t0*t2*t51*t8 + (803.0/64.0)*t0*t4*t51*t8 + (45.0/4.0)*t0*t51*t7*t8 - 105.0/2.0*t1*t11 + (149.0/16.0)*t1*t12*t7 - 135.0/4.0*t1*t13 + (21.0/2.0)*t1*t14*t15 + (9.0/4.0)*t1*t15*t17 - 45.0/8.0*t1*t16 - t1*t21 + (475.0/64.0)*t1*t4*t8 + (447.0/16.0)*t1*t4 + (405.0/4.0)*t1*t7*t8 + (225.0/4.0)*t1*t7*t9 + (135.0/256.0)*t1 + 45*t10*t12*t22*t43 + 105*t10*t12*t23*t24 + (15.0/2.0)*t10*t14*t23*t24 + (33.0/4.0)*t10*t14*t29 + 210*t10*t22*t43*t9 + 210*t10*t23*t24*t9 - 315*t11*t29 - 105.0/2.0*t11*t52 + (87.0/4.0)*t12*t22*t43*t7 - t13*t37 - t13*t53 + 63*t14*t15*t29 + (15.0/2.0)*t15*t17*t29 - 33.0/2.0*t16*t44 + 2*t18*t19*t22*t43 + 2*t18*t19*t23*t24 + (335.0/8.0)*t2*t22*t43*t9 + 45*t2*t22*t43 + (15.0/2.0)*t2*t23*t24*t8 + 105*t2*t23*t24 + (639.0/4.0)*t2*t29*t8 - t2*t37 - t2*t53 - 3*t20*t29 - t21*t52 + (261.0/4.0)*t22*t4*t43 + (693.0/16.0)*t22*t43*t5 + 165*t22*t43*t7*t9 + (1377.0/64.0)*t22*t43 + (15.0/16.0)*t23*t24*t4*t8 - t25*t27 - 315*t25*t30 - 75*t25*t31 - 17.0/8.0*t25*t32 - 5.0/4.0*t25*t33 - t25*t35 - 7*t25*t36 - 15.0/4.0*t25*t4 - 15.0/16.0*t25*t5 - 3.0/64.0*t25 - t27*t44 - 165.0/2.0*t29*t31 - 81.0/16.0*t29*t32 - 213.0/8.0*t29*t33 + (789.0/32.0)*t29*t4*t8 + (675.0/2.0)*t29*t7*t8 - 783.0/128.0*t29 - 447.0/8.0*t3*t8 - 323.0/32.0*t3*t9 - 135.0/4.0*t3 - 135*t30*t44 - 255.0/4.0*t31*t52 - 803.0/32.0*t32*t52 - 51.0/16.0*t33*t52 - t35*t44 - 639.0/8.0*t38 - 789.0/32.0*t39 - 153.0/16.0*t4*t52 - t42*t50 - 495.0/2.0*t44*t45 - 1005.0/8.0*t44*t46 - 1071.0/32.0*t44*t47 - 261.0/2.0*t44*t48 - 693.0/16.0*t44*t49 - 803.0/64.0*t5*t52 - 3825.0/256.0*t52 - 475.0/64.0*t6);
    o_mean_mean_mean_var[i] = (1.0/512.0)*t54*(10080*C*t2*t24*t56 + 10080*C*t2*t29*t60 + 24480*C*t2*t51*t55 + 40320*C*t29*t60*t7 + 19272*C*t4*t51*t55 + 7650*C*t5*t51*t55 + 13440*C*t51*t55*t7 + 8640*m*t10*t12*t43 + 40320*m*t10*t43*t9 + 4080*m*t12*t43*t7 + 384*m*t18*t19*t43 + 6792*m*t2*t43*t9 + 8640*m*t2*t43 + 12240*m*t4*t43 + 7692*m*t43*t5 + 29760*m*t43*t7*t9 + 3213*m*t43 + 6720*t10*t12*t24*t56 + 672*t10*t14*t24*t56 + 672*t10*t14*t29*t60 + 1632*t10*t14*t51*t55 + 13440*t10*t24*t56*t9 - 40320*t11*t61 - t13*t62 + 8064*t14*t15*t29*t60 + 2688*t14*t15*t51*t55 + 960*t15*t17*t29*t60 + 64*t15*t17*t51*t55 - 2976*t16*t65 + 128*t18*t19*t24*t56 + 3360*t2*t24*t56*t8 + 152*t2*t24*t56*t9 + 6720*t2*t24*t56 + 16992*t2*t29*t60*t8 + 504*t2*t29*t60*t9 + 4896*t2*t51*t55*t8 - t2*t62 - 384*t20*t61 + 108*t24*t4*t56*t8 + 15*t24*t56 - 40320*t26*t65 + 1620*t29*t4*t60*t8 + 43200*t29*t60*t7*t8 - 20160*t30*t57 - 25920*t30*t65 - 6720*t31*t57 - 6720*t31*t61 - 560*t33*t57 - 2832*t33*t61 - 8064*t34*t65 - 448*t36*t57 - 576*t36*t65 - 8496*t38*t60 - 1620*t39*t60 + 3212*t4*t51*t55*t8 - 1680*t4*t57 - 44640*t45*t65 - 456*t46*t57 - 20376*t46*t65 - 82*t47*t57 - 4590*t47*t65 - 24480*t48*t65 - 7692*t49*t65 - 108*t5*t57 + 2880*t51*t55*t7*t8 - t55*t66 - t57*t58 - t57*t59 - 1512*t60*t63 - 162*t60*t64 - 243*t61 - t67*t68 - t68*t69 - t68*t70 - t68*t71 - t68*t72 - t68*t73 - t68*t74 - t68*t75 - t68*t76);
    o_mean_mean_var_var[i] = (1.0/1024.0)*t54*(24480*C*t2*t22*t51 + 1440*C*t2*t23*t29 + 19272*C*t22*t4*t51 + 7650*C*t22*t5*t51 + 13440*C*t22*t51*t7 + 13440*C*t23*t29*t7 + 5760*t10*t12*t43 + 1632*t10*t14*t22*t51 + 96*t10*t14*t23*t29 + 26880*t10*t43*t9 + 2656*t12*t43*t7 - t13*t78 + 2688*t14*t15*t22*t51 + 2688*t14*t15*t23*t29 + 64*t15*t17*t22*t51 + 320*t15*t17*t23*t29 - 1856*t16*t43 + 256*t18*t19*t43 + 4896*t2*t22*t51*t8 + 4512*t2*t23*t29*t8 + 488*t2*t23*t29*t9 + 3696*t2*t43*t9 + 5760*t2*t43 - t2*t78 + 3212*t22*t4*t51*t8 + 2880*t22*t51*t7*t8 - t22*t66 + 60*t23*t29*t4*t8 + 14400*t23*t29*t7*t8 + 27*t23*t29 - 2256*t23*t38 - 60*t23*t39 - 1464*t23*t63 - 126*t23*t64 - 26880*t26*t43 - 17280*t30*t43 - 960*t31*t77 - 752*t33*t77 - 5376*t34*t43 + 7968*t4*t43 - 1836*t43*t47 - 4712*t43*t49 + 4712*t43*t5 + 18560*t43*t7*t9 + 1530*t43 - 384*t50 - t67*t84 - t69*t84 - t70*t84 - t71*t77 - t71*t84 - t72*t84 - t73*t84 - t74*t84 - t75*t84 - t76*t77 - t76*t84 - 27840*t80 - 11088*t82 - 15936*t83);
    o_mean_var_var_var[i] = (1.0/2048.0)*t54*(24480*C*m*t2*t51 + 19272*C*m*t4*t51 + 7650*C*m*t5*t51 + 13440*C*m*t51*t7 + 1632*m*t10*t14*t51 + 2688*m*t14*t15*t51 + 64*m*t15*t17*t51 + 4896*m*t2*t51*t8 + 3212*m*t4*t51*t8 + 2880*m*t51*t7*t8 - m*t66 + 2880*t10*t12*t43*t60 + 13440*t10*t43*t60*t9 + 1296*t12*t43*t60*t7 - 864*t16*t85 + 128*t18*t19*t43*t60 + 1432*t2*t43*t60*t9 + 2880*t2*t43*t60 - 8640*t30*t85 - 192*t36*t85 + 3888*t4*t43*t60 + 2148*t43*t5*t60 + 8640*t43*t60*t7*t9 + 459*t43*t60 - 306*t47*t85 - t58*t85 - t59*t85 - 2148*t60*t8*t81 - 12960*t60*t80 - 4296*t60*t82 - 7776*t60*t83 - t67*t86 - t69*t86 - t70*t86 - t71*t86 - t72*t86 - t73*t86 - t74*t86 - t75*t86 - t76*t86);
    o_var_var_var_var[i] = (1.0/4096.0)*t51*t54*(24480*C*t2 + 19272*C*t4 + 7650*C*t5 + 13440*C*t7 + 1632*t10*t14 + 2688*t14*t15 + 64*t15*t17 + 4896*t2*t8 + 3212*t4*t8 - t67 - t69 + 2880*t7*t8 - t70 - t71 - t72 - t73 - t74 - t75 - t76 - 3825);
  }
  (void) S;
  return List::create(Named("mean_mean_mean_mean") = o_mean_mean_mean_mean, Named("mean_mean_mean_var") = o_mean_mean_mean_var, Named("mean_mean_var_var") = o_mean_mean_var_var, Named("mean_var_var_var") = o_mean_var_var_var, Named("var_var_var_var") = o_var_var_var_var);
}

// 1 in y and 1 in (mean, var)
// [[Rcpp::export]]
List lognormal2_cross_y_cpp(NumericVector y, NumericVector mean, NumericVector var) {
  const int n = y.size();
  const int n_mean = mean.size(), n_var = var.size();
  NumericVector o_mean(n);
  NumericVector o_var(n);
  double m = 0.0, v = 0.0, S = 0.0;
  double t0 = 0.0;
  double t2 = 0.0;
  double t4 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; v = var[i % n_var];
    S = std::log1p(v / (m * m));
    t0 = 1.0/S;
    t2 = 1.0/m;
    t4 = 1.0/(std::pow(m, 2) + v);
  };
  const bool scalar = n_mean == 1 && n_var == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t1 = 1.0/yy;
    const double t3 = 2*C*t0;
    o_mean[i] = t0*t1*(m*t3*t4 - t2*t3 + t2);
    o_var[i] = C*t1*t4/std::pow(S, 2);
  }
  (void) S;
  return List::create(Named("mean") = o_mean, Named("var") = o_var);
}

// 2 in y and 1 in (mean, var)
// [[Rcpp::export]]
List lognormal2_cross2_y_cpp(NumericVector y, NumericVector mean, NumericVector var) {
  const int n = y.size();
  const int n_mean = mean.size(), n_var = var.size();
  NumericVector o_mean(n);
  NumericVector o_var(n);
  double m = 0.0, v = 0.0, S = 0.0;
  double t0 = 0.0;
  double t2 = 0.0;
  double t3 = 0.0;
  double t4 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; v = var[i % n_var];
    S = std::log1p(v / (m * m));
    t0 = 1.0/S;
    t2 = 1.0/m;
    t3 = 2*t0;
    t4 = 1.0/(std::pow(m, 2) + v);
  };
  const bool scalar = n_mean == 1 && n_var == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t1 = std::pow(yy, -2);
    o_mean[i] = t0*t1*(-C*m*t3*t4 + 2*C*t0*t2 + 2*m*t0*t4 - t2*t3 - t2);
    o_var[i] = t1*t4*(1 - C)/std::pow(S, 2);
  }
  (void) S;
  return List::create(Named("mean") = o_mean, Named("var") = o_var);
}

// 1 in y and 2 in (mean, var)
// [[Rcpp::export]]
List lognormal2_grad_y_hess_cpp(NumericVector y, NumericVector mean, NumericVector var) {
  const int n = y.size();
  const int n_mean = mean.size(), n_var = var.size();
  NumericVector o_mean_mean(n);
  NumericVector o_var_var(n);
  NumericVector o_mean_var(n);
  double m = 0.0, v = 0.0, S = 0.0;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double t3 = 0.0;
  double t4 = 0.0;
  double t5 = 0.0;
  double t6 = 0.0;
  double t11 = 0.0;
  double t13 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; v = var[i % n_var];
    S = std::log1p(v / (m * m));
    t0 = std::pow(m, 2);
    t1 = 1.0/t0;
    t2 = 1.0/S;
    t3 = std::pow(S, -2);
    t4 = 8*t3;
    t5 = 1.0/(t0 + v);
    t6 = 4*t2;
    t11 = 2*t2;
    t13 = 1.0/m;
  };
  const bool scalar = n_mean == 1 && n_var == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t7 = C/(std::pow(m, 4) + 2*t0*v + std::pow(v, 2));
    const double t8 = t0*t7;
    const double t9 = 1.0/yy;
    const double t10 = t2*t9;
    const double t12 = t3*t7;
    o_mean_mean[i] = t10*(2*C*t1*t2 - C*t1*t4 + 2*C*t2*t5 + 16*C*t3*t5 + 4*t1*t2 - t1 - t4*t8 - t5*t6 - t6*t8);
    o_var_var[i] = -t12*t9*(t11 + 1);
    o_mean_var[i] = t10*(4*C*t13*t3*t5 - m*t11*t7 - 4*m*t12 - t13*t2*t5);
  }
  (void) S;
  return List::create(Named("mean_mean") = o_mean_mean, Named("var_var") = o_var_var, Named("mean_var") = o_mean_var);
}

// 2 in y and 2 in (mean, var)
// [[Rcpp::export]]
List lognormal2_hess_y_hess_cpp(NumericVector y, NumericVector mean, NumericVector var) {
  const int n = y.size();
  const int n_mean = mean.size(), n_var = var.size();
  NumericVector o_mean_mean(n);
  NumericVector o_var_var(n);
  NumericVector o_mean_var(n);
  double m = 0.0, v = 0.0, S = 0.0;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double t3 = 0.0;
  double t4 = 0.0;
  double t5 = 0.0;
  double t6 = 0.0;
  double t7 = 0.0;
  double t8 = 0.0;
  double t10 = 0.0;
  double t11 = 0.0;
  double t14 = 0.0;
  double t15 = 0.0;
  double t16 = 0.0;
  double t17 = 0.0;
  double t18 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; v = var[i % n_var];
    S = std::log1p(v / (m * m));
    t0 = std::pow(m, 2);
    t1 = 1.0/t0;
    t2 = std::pow(S, -2);
    t3 = 8*t2;
    t4 = 1.0/S;
    t5 = 2*t4;
    t6 = t1*t5;
    t7 = 1.0/(t0 + v);
    t8 = t2*t7;
    t10 = 1.0/(std::pow(m, 4) + 2*t0*v + std::pow(v, 2));
    t11 = t0*t10;
    t14 = t10*t2;
    t15 = 1.0/m;
    t16 = 4*t15*t8;
    t17 = 4*m*t14;
    t18 = m*t10;
  };
  const bool scalar = n_mean == 1 && n_var == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t9 = C*t5;
    const double t12 = std::pow(yy, -2);
    const double t13 = t12*t4;
    o_mean_mean[i] = t13*(8*C*t0*t10*t2 + 4*C*t0*t10*t4 + 8*C*t1*t2 - C*t6 - 16*C*t8 - t1*t3 + t1 - t11*t3 - 4*t11*t4 + 16*t2*t7 + 6*t4*t7 - t6 - t7*t9);
    o_var_var[i] = t12*t14*(C - t5 + t9 - 1);
    o_mean_var[i] = t13*(-C*t16 + C*t17 + t15*t4*t7 + t16 - t17 - t18*t5 + t18*t9);
  }
  (void) S;
  return List::create(Named("mean_mean") = o_mean_mean, Named("var_var") = o_var_var, Named("mean_var") = o_mean_var);
}

// order 1 of log f in y
// [[Rcpp::export]]
List lognormal2_dy1_cpp(NumericVector y, NumericVector mean, NumericVector var) {
  const int n = y.size();
  const int n_mean = mean.size(), n_var = var.size();
  NumericVector o_y(n);
  double m = 0.0, v = 0.0, S = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; v = var[i % n_var];
    S = std::log1p(v / (m * m));
  };
  const bool scalar = n_mean == 1 && n_var == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    o_y[i] = -(C/S + 3.0/2.0)/yy;
  }
  (void) S;
  return List::create(Named("y") = o_y);
}

// order 2 of log f in y
// [[Rcpp::export]]
List lognormal2_dy2_cpp(NumericVector y, NumericVector mean, NumericVector var) {
  const int n = y.size();
  const int n_mean = mean.size(), n_var = var.size();
  NumericVector o_y(n);
  double m = 0.0, v = 0.0, S = 0.0;
  double t0 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; v = var[i % n_var];
    S = std::log1p(v / (m * m));
    t0 = 1.0/S;
  };
  const bool scalar = n_mean == 1 && n_var == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    o_y[i] = (C*t0 - t0 + 3.0/2.0)/std::pow(yy, 2);
  }
  (void) S;
  return List::create(Named("y") = o_y);
}

// order 3 of log f in y
// [[Rcpp::export]]
List lognormal2_dy3_cpp(NumericVector y, NumericVector mean, NumericVector var) {
  const int n = y.size();
  const int n_mean = mean.size(), n_var = var.size();
  NumericVector o_y(n);
  double m = 0.0, v = 0.0, S = 0.0;
  double t0 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; v = var[i % n_var];
    S = std::log1p(v / (m * m));
    t0 = 1.0/S;
  };
  const bool scalar = n_mean == 1 && n_var == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    o_y[i] = (-2*C*t0 + 3*t0 - 3)/std::pow(yy, 3);
  }
  (void) S;
  return List::create(Named("y") = o_y);
}

// order 4 of log f in y
// [[Rcpp::export]]
List lognormal2_dy4_cpp(NumericVector y, NumericVector mean, NumericVector var) {
  const int n = y.size();
  const int n_mean = mean.size(), n_var = var.size();
  NumericVector o_y(n);
  double m = 0.0, v = 0.0, S = 0.0;
  double t0 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; v = var[i % n_var];
    S = std::log1p(v / (m * m));
    t0 = 1.0/S;
  };
  const bool scalar = n_mean == 1 && n_var == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    o_y[i] = (6*C*t0 - 11*t0 + 9)/std::pow(yy, 4);
  }
  (void) S;
  return List::create(Named("y") = o_y);
}
