#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// The Weibull in its mean m and shape s. With b = m / Gamma(1 + 1/s),
// C = log(y/m) and w = C + lgamma(1 + 1/s) = log(y/b), the log-density is
//   log s - log y + s w - exp(s w).
// Every component below is a closed form derived offline with sympy
// (stabilita/gen_weibull3.py), one exported function per order and per
// surface; P_k is psigamma(1 + 1/s, k) and LG is lgamma(1 + 1/s), both
// computed once per parameter value. Each is checked numerically in the
// tests.


// order 1 of log f in (mean, sigma)
// [[Rcpp::export]]
List weibull3_gradient_cpp(NumericVector y, NumericVector mean, NumericVector sigma) {
  const int n = y.size();
  const int n_mean = mean.size(), n_sigma = sigma.size();
  NumericVector o_mean(n);
  NumericVector o_sigma(n);
  double m = 0.0, s = 0.0, LG = 0.0;
  double P0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; s = sigma[i % n_sigma];
    const double x1 = 1.0 + 1.0 / s; (void) x1;
    LG = R::lgammafn(x1);
    P0 = R::psigamma(x1, 0.0);
    t1 = 1.0/s;
    t2 = P0*t1;
  };
  const bool scalar = n_mean == 1 && n_sigma == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t0 = std::exp(C*s)*std::exp(LG*s);
    o_mean[i] = s*(t0 - 1)/m;
    o_sigma[i] = -C*t0 + C - LG*t0 + LG + t0*t2 + t1 - t2;
  }
  (void) LG;
  return List::create(Named("mean") = o_mean, Named("sigma") = o_sigma);
}

// order 2 of log f in (mean, sigma)
// [[Rcpp::export]]
List weibull3_hessian_cpp(NumericVector y, NumericVector mean, NumericVector sigma) {
  const int n = y.size();
  const int n_mean = mean.size(), n_sigma = sigma.size();
  NumericVector o_mean_mean(n);
  NumericVector o_sigma_sigma(n);
  NumericVector o_mean_sigma(n);
  double m = 0.0, s = 0.0, LG = 0.0;
  double P0 = 0.0;
  double P1 = 0.0;
  double t2 = 0.0;
  double t3 = 0.0;
  double t6 = 0.0;
  double t7 = 0.0;
  double t8 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; s = sigma[i % n_sigma];
    const double x1 = 1.0 + 1.0 / s; (void) x1;
    LG = R::lgammafn(x1);
    P0 = R::psigamma(x1, 0.0);
    P1 = R::psigamma(x1, 1.0);
    t2 = LG*s;
    t3 = std::exp(t2);
    t6 = std::pow(s, -2);
    t7 = std::pow(s, -3);
    t8 = 1.0/s;
  };
  const bool scalar = n_mean == 1 && n_sigma == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t0 = C*s;
    const double t1 = std::exp(t0);
    const double t4 = t1*t3;
    const double t5 = t4 - 1;
    o_mean_mean[i] = s*(-s*t4 - t5)/std::pow(m, 2);
    o_sigma_sigma[i] = -std::pow(C, 2)*t4 - 2*C*LG*t4 + 2*C*P0*t1*t3*t8 - std::pow(LG, 2)*t4 + 2*LG*P0*t1*t3*t8 - std::pow(P0, 2)*t4*t6 - P1*t4*t7 + P1*t7 - t6;
    o_mean_sigma[i] = (-P0*t4 + t0*t4 + t2*t4 + t5)/m;
  }
  (void) LG;
  return List::create(Named("mean_mean") = o_mean_mean, Named("sigma_sigma") = o_sigma_sigma, Named("mean_sigma") = o_mean_sigma);
}

// order 2, expected
// [[Rcpp::export]]
List weibull3_expected_hessian_cpp(NumericVector y, NumericVector mean, NumericVector sigma) {
  const int n = y.size();
  const int n_mean = mean.size(), n_sigma = sigma.size();
  NumericVector o_mean_mean(n);
  NumericVector o_sigma_sigma(n);
  NumericVector o_mean_sigma(n);
  double m = 0.0, s = 0.0, LG = 0.0;
  double P0 = 0.0;
  double t0 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; s = sigma[i % n_sigma];
    const double x1 = 1.0 + 1.0 / s; (void) x1;
    LG = R::lgammafn(x1);
    P0 = R::psigamma(x1, 0.0);
    t0 = std::pow(s, 2);
  };
  const bool scalar = n_mean == 1 && n_sigma == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    o_mean_mean[i] = -1.0*t0/std::pow(m, 2);
    o_sigma_sigma[i] = (-1.0*std::pow(P0, 2) + 0.84556867019693428*P0 - 1.8236806608528794)/t0;
    o_mean_sigma[i] = (0.42278433509846714 - 1.0*P0)/m;
  }
  (void) LG;
  return List::create(Named("mean_mean") = o_mean_mean, Named("sigma_sigma") = o_sigma_sigma, Named("mean_sigma") = o_mean_sigma);
}

// order 3 of log f in (mean, sigma)
// [[Rcpp::export]]
List weibull3_deriv3_cpp(NumericVector y, NumericVector mean, NumericVector sigma) {
  const int n = y.size();
  const int n_mean = mean.size(), n_sigma = sigma.size();
  NumericVector o_mean_mean_mean(n);
  NumericVector o_mean_mean_sigma(n);
  NumericVector o_mean_sigma_sigma(n);
  NumericVector o_sigma_sigma_sigma(n);
  double m = 0.0, s = 0.0, LG = 0.0;
  double P0 = 0.0;
  double P1 = 0.0;
  double P2 = 0.0;
  double t2 = 0.0;
  double t3 = 0.0;
  double t7 = 0.0;
  double t10 = 0.0;
  double t12 = 0.0;
  double t13 = 0.0;
  double t14 = 0.0;
  double t15 = 0.0;
  double t16 = 0.0;
  double t17 = 0.0;
  double t18 = 0.0;
  double t19 = 0.0;
  double t23 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; s = sigma[i % n_sigma];
    const double x1 = 1.0 + 1.0 / s; (void) x1;
    LG = R::lgammafn(x1);
    P0 = R::psigamma(x1, 0.0);
    P1 = R::psigamma(x1, 1.0);
    P2 = R::psigamma(x1, 2.0);
    t2 = LG*s;
    t3 = std::exp(t2);
    t7 = std::pow(s, 2);
    t10 = 2*LG;
    t12 = std::pow(LG, 2);
    t13 = 1.0/t7;
    t14 = 1.0/s;
    t15 = std::pow(P0, 2);
    t16 = std::pow(s, -3);
    t17 = std::pow(s, -5);
    t18 = std::pow(s, -4);
    t19 = 3*P1;
    t23 = t13*t15;
  };
  const bool scalar = n_mean == 1 && n_sigma == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t0 = C*s;
    const double t1 = std::exp(t0);
    const double t4 = t1*t3;
    const double t5 = 2*t4;
    const double t6 = 3*t4;
    const double t8 = t4*t7;
    const double t9 = 2*C;
    const double t11 = std::pow(C, 2);
    const double t20 = C*t6;
    const double t21 = LG*t6;
    const double t22 = t16*t19*t4;
    o_mean_mean_mean[i] = s*(s*t6 + t5 + t8 - 2)/std::pow(m, 3);
    o_mean_mean_sigma[i] = (-C*t8 - LG*t8 + P0*s*t1*t3 + P0*t1*t3 - s*t5 - t0*t4 - t2*t4 - t4 + 1)/std::pow(m, 2);
    o_mean_sigma_sigma[i] = t4*(-P0*t10 - 2*P0*t14 - P0*t9 + P1*t13 + s*t11 + s*t12 + t0*t10 + t10 + t14*t15 + t9)/m;
    o_sigma_sigma_sigma[i] = -std::pow(C, 3)*t4 + 6*C*LG*P0*t1*t14*t3 - C*t22 - std::pow(LG, 3)*t4 - LG*t22 + std::pow(P0, 3)*t1*t16*t3 + 3*P0*P1*t1*t18*t3 + 3*P0*t1*t11*t14*t3 + 3*P0*t1*t12*t14*t3 + 3*P1*t1*t18*t3 + P2*t1*t17*t3 - P2*t17 - t11*t21 - t12*t20 + 2*t16 - t18*t19 - t20*t23 - t21*t23;
  }
  (void) LG;
  return List::create(Named("mean_mean_mean") = o_mean_mean_mean, Named("mean_mean_sigma") = o_mean_mean_sigma, Named("mean_sigma_sigma") = o_mean_sigma_sigma, Named("sigma_sigma_sigma") = o_sigma_sigma_sigma);
}

// order 3, expected
// [[Rcpp::export]]
List weibull3_deriv3_expected_cpp(NumericVector y, NumericVector mean, NumericVector sigma) {
  const int n = y.size();
  const int n_mean = mean.size(), n_sigma = sigma.size();
  NumericVector o_mean_mean_mean(n);
  NumericVector o_mean_mean_sigma(n);
  NumericVector o_mean_sigma_sigma(n);
  NumericVector o_sigma_sigma_sigma(n);
  double m = 0.0, s = 0.0, LG = 0.0;
  double P0 = 0.0;
  double P1 = 0.0;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double t3 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; s = sigma[i % n_sigma];
    const double x1 = 1.0 + 1.0 / s; (void) x1;
    LG = R::lgammafn(x1);
    P0 = R::psigamma(x1, 0.0);
    P1 = R::psigamma(x1, 1.0);
    t0 = 1.0*s;
    t1 = 1.0/s;
    t2 = std::pow(P0, 2);
    t3 = P1*t1;
  };
  const bool scalar = n_mean == 1 && n_sigma == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    o_mean_mean_mean[i] = std::pow(s, 2)*(t0 + 3.0)/std::pow(m, 3);
    o_mean_mean_sigma[i] = (P0*t0 + 1.0*P0 - 2.4227843350984671*s - 0.42278433509846714)/std::pow(m, 2);
    o_mean_sigma_sigma[i] = t1*(-2.8455686701969343*P0 + 1.0*t2 + 1.0*t3 + 1.6692493310498137)/m;
    o_sigma_sigma_sigma[i] = (1.0*std::pow(P0, 3) + 3.0*P0*t3 + 2.4710419825586382*P0 - 1.2683530052954014*t2 - 1.2683530052954014*t3 + 1.5105384845174824)/std::pow(s, 3);
  }
  (void) LG;
  return List::create(Named("mean_mean_mean") = o_mean_mean_mean, Named("mean_mean_sigma") = o_mean_mean_sigma, Named("mean_sigma_sigma") = o_mean_sigma_sigma, Named("sigma_sigma_sigma") = o_sigma_sigma_sigma);
}

// order 4 of log f in (mean, sigma)
// [[Rcpp::export]]
List weibull3_deriv4_cpp(NumericVector y, NumericVector mean, NumericVector sigma) {
  const int n = y.size();
  const int n_mean = mean.size(), n_sigma = sigma.size();
  NumericVector o_mean_mean_mean_mean(n);
  NumericVector o_mean_mean_mean_sigma(n);
  NumericVector o_mean_mean_sigma_sigma(n);
  NumericVector o_mean_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_sigma_sigma(n);
  double m = 0.0, s = 0.0, LG = 0.0;
  double P0 = 0.0;
  double P1 = 0.0;
  double P2 = 0.0;
  double P3 = 0.0;
  double t2 = 0.0;
  double t3 = 0.0;
  double t7 = 0.0;
  double t9 = 0.0;
  double t11 = 0.0;
  double t15 = 0.0;
  double t16 = 0.0;
  double t18 = 0.0;
  double t19 = 0.0;
  double t20 = 0.0;
  double t21 = 0.0;
  double t22 = 0.0;
  double t26 = 0.0;
  double t28 = 0.0;
  double t29 = 0.0;
  double t30 = 0.0;
  double t32 = 0.0;
  double t33 = 0.0;
  double t34 = 0.0;
  double t35 = 0.0;
  double t36 = 0.0;
  double t37 = 0.0;
  double t38 = 0.0;
  double t39 = 0.0;
  double t41 = 0.0;
  double t43 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; s = sigma[i % n_sigma];
    const double x1 = 1.0 + 1.0 / s; (void) x1;
    LG = R::lgammafn(x1);
    P0 = R::psigamma(x1, 0.0);
    P1 = R::psigamma(x1, 1.0);
    P2 = R::psigamma(x1, 2.0);
    P3 = R::psigamma(x1, 3.0);
    t2 = LG*s;
    t3 = std::exp(t2);
    t7 = std::pow(s, 3);
    t9 = std::pow(s, 2);
    t11 = 3*P0;
    t15 = 2*LG;
    t16 = std::pow(P0, 2);
    t18 = std::pow(LG, 2);
    t19 = 1.0/t9;
    t20 = P1*t19;
    t21 = 1.0/s;
    t22 = t16*t21;
    t26 = 3*t18;
    t28 = std::pow(LG, 3);
    t29 = std::pow(s, -4);
    t30 = P0*t21;
    t32 = 3*LG;
    t33 = 1.0/t7;
    t34 = P1*t33;
    t35 = t16*t19;
    t36 = std::pow(P0, 3);
    t37 = std::pow(s, -5);
    t38 = std::pow(s, -6);
    t39 = std::pow(s, -7);
    t41 = P1*t37;
    t43 = P2*t38;
  };
  const bool scalar = n_mean == 1 && n_sigma == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t0 = C*s;
    const double t1 = std::exp(t0);
    const double t4 = t1*t3;
    const double t5 = 6*t4;
    const double t6 = s*t4;
    const double t8 = t4*t7;
    const double t10 = 2*t4;
    const double t12 = t4*t9;
    const double t13 = 3*t12;
    const double t14 = 2*C;
    const double t17 = std::pow(C, 2);
    const double t23 = 6*C;
    const double t24 = LG*t23;
    const double t25 = 3*t17;
    const double t27 = std::pow(C, 3);
    const double t31 = 3*C;
    const double t40 = 4*t4;
    const double t42 = 12*t4*t41;
    const double t44 = t17*t5;
    const double t45 = 12*C*LG*t4;
    const double t46 = t18*t5;
    o_mean_mean_mean_mean[i] = s*(-t5*t9 - t5 - 11*t6 - t8 + 6)/std::pow(m, 4);
    o_mean_mean_mean_sigma[i] = (C*t13 + C*t8 + LG*t13 + LG*t8 - P0*t10 - P0*t12 + s*t5 + t0*t10 + t10*t2 + t10 - t11*t6 + t13 - 2)/std::pow(m, 3);
    o_mean_mean_sigma_sigma[i] = t4*(2*C*P0*s + 2*C*P0 + 2*LG*P0*s + 2*LG*P0 - LG*t14*t9 + 2*P0*t21 + 4*P0 - P1*t21 - s*t17 - s*t18 - t0*t15 - 4*t0 - t14 - t15 - t16 - t17*t9 - t18*t9 - 4*t2 - t20 - t22 - 2)/std::pow(m, 2);
    o_mean_sigma_sigma_sigma[i] = t4*(-6*LG*t30 - P0*t24 - P0*t25 - P0*t26 - P2*t29 + s*t27 + s*t28 + t0*t26 - t11*t34 - t19*t36 + t2*t25 + t20*t31 + t20*t32 + t22*t31 + t22*t32 - t23*t30 + t24 + t25 + t26 + 3*t35)/m;
    o_sigma_sigma_sigma_sigma[i] = -std::pow(C, 4)*t4 + 12*C*P0*P1*t1*t29*t3 + 12*C*P0*t1*t18*t21*t3 + 12*C*P1*t1*t29*t3 + 4*C*P2*t1*t3*t37 + 4*C*t1*t3*t33*t36 - C*t28*t40 - std::pow(LG, 4)*t4 + 12*LG*P0*P1*t1*t29*t3 + 12*LG*P0*t1*t17*t21*t3 + 12*LG*P1*t1*t29*t3 + 4*LG*P2*t1*t3*t37 + 4*LG*t1*t3*t33*t36 - LG*t27*t40 - std::pow(P0, 4)*t29*t4 + 4*P0*t1*t21*t27*t3 + 4*P0*t1*t21*t28*t3 - 4*P0*t4*t43 - P0*t42 - 3*std::pow(P1, 2)*t38*t4 + 12*P1*t37 + 8*P2*t38 - P3*t39*t4 + P3*t39 - t16*t41*t5 - t18*t44 - 6*t29 - t34*t44 - t34*t45 - t34*t46 - t35*t44 - t35*t45 - t35*t46 - 8*t4*t43 - t42;
  }
  (void) LG;
  return List::create(Named("mean_mean_mean_mean") = o_mean_mean_mean_mean, Named("mean_mean_mean_sigma") = o_mean_mean_mean_sigma, Named("mean_mean_sigma_sigma") = o_mean_mean_sigma_sigma, Named("mean_sigma_sigma_sigma") = o_mean_sigma_sigma_sigma, Named("sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma);
}

// order 4, expected
// [[Rcpp::export]]
List weibull3_deriv4_expected_cpp(NumericVector y, NumericVector mean, NumericVector sigma) {
  const int n = y.size();
  const int n_mean = mean.size(), n_sigma = sigma.size();
  NumericVector o_mean_mean_mean_mean(n);
  NumericVector o_mean_mean_mean_sigma(n);
  NumericVector o_mean_mean_sigma_sigma(n);
  NumericVector o_mean_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_sigma_sigma(n);
  double m = 0.0, s = 0.0, LG = 0.0;
  double P0 = 0.0;
  double P1 = 0.0;
  double P2 = 0.0;
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
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; s = sigma[i % n_sigma];
    const double x1 = 1.0 + 1.0 / s; (void) x1;
    LG = R::lgammafn(x1);
    P0 = R::psigamma(x1, 0.0);
    P1 = R::psigamma(x1, 1.0);
    P2 = R::psigamma(x1, 2.0);
    t0 = std::pow(s, 2);
    t1 = 1.0*t0;
    t2 = 3.0*P0;
    t3 = std::pow(P0, 2);
    t4 = 1.0*t3;
    t5 = 1.0/s;
    t6 = 1.0/t0;
    t7 = 1.0*P1;
    t8 = std::pow(P0, 3);
    t9 = P2*t6;
    t10 = P1*t5;
  };
  const bool scalar = n_mean == 1 && n_sigma == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    o_mean_mean_mean_mean[i] = -t0*(6.0*s + t1 + 11.0)/std::pow(m, 4);
    o_mean_mean_mean_sigma[i] = (-P0*t1 - 2.0*P0 - s*t2 + 7.2683530052954014*s + 3.4227843350984671*t0 + 0.84556867019693428)/std::pow(m, 3);
    o_mean_mean_sigma_sigma[i] = (2.8455686701969343*P0*t5 + 4.8455686701969343*P0 - t4*t5 - t4 - t5*t7 - 1.6692493310498137*t5 - t6*t7 - 4.5148180012467479)/std::pow(m, 2);
    o_mean_sigma_sigma_sigma[i] = t6*(-5.007747993149441*P0 + 1.2683530052954014*P1*t5 - t10*t2 + 4.2683530052954014*t3 - 1.0*t8 - 1.0*t9 + 2.9605034980411558)/m;
    o_sigma_sigma_sigma_sigma[i] = (-1.0*std::pow(P0, 4) - 6.9265879788183943*P0*P1*t5 - 4.0*P0*t9 + 1.9578460619300704*P0 - 3.0*std::pow(P1, 2)*t6 + 0.13132805606432934*P1*t5 + 1.6911373403938686*P2*t6 - 6.0*t10*t3 - 4.9420839651172763*t3 + 1.6911373403938686*t8 - 7.7819762580843336)/std::pow(s, 4);
  }
  (void) LG;
  return List::create(Named("mean_mean_mean_mean") = o_mean_mean_mean_mean, Named("mean_mean_mean_sigma") = o_mean_mean_mean_sigma, Named("mean_mean_sigma_sigma") = o_mean_mean_sigma_sigma, Named("mean_sigma_sigma_sigma") = o_mean_sigma_sigma_sigma, Named("sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma);
}

// order 5 of log f in (mean, sigma)
// [[Rcpp::export]]
List weibull3_deriv5_cpp(NumericVector y, NumericVector mean, NumericVector sigma) {
  const int n = y.size();
  const int n_mean = mean.size(), n_sigma = sigma.size();
  NumericVector o_mean_mean_mean_mean_mean(n);
  NumericVector o_mean_mean_mean_mean_sigma(n);
  NumericVector o_mean_mean_mean_sigma_sigma(n);
  NumericVector o_mean_mean_sigma_sigma_sigma(n);
  NumericVector o_mean_sigma_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_sigma_sigma_sigma(n);
  double m = 0.0, s = 0.0, LG = 0.0;
  double P0 = 0.0;
  double P1 = 0.0;
  double P2 = 0.0;
  double P3 = 0.0;
  double P4 = 0.0;
  double t2 = 0.0;
  double t3 = 0.0;
  double t6 = 0.0;
  double t8 = 0.0;
  double t10 = 0.0;
  double t15 = 0.0;
  double t17 = 0.0;
  double t18 = 0.0;
  double t19 = 0.0;
  double t21 = 0.0;
  double t22 = 0.0;
  double t23 = 0.0;
  double t26 = 0.0;
  double t27 = 0.0;
  double t28 = 0.0;
  double t29 = 0.0;
  double t30 = 0.0;
  double t31 = 0.0;
  double t32 = 0.0;
  double t34 = 0.0;
  double t37 = 0.0;
  double t38 = 0.0;
  double t39 = 0.0;
  double t41 = 0.0;
  double t44 = 0.0;
  double t45 = 0.0;
  double t46 = 0.0;
  double t47 = 0.0;
  double t48 = 0.0;
  double t49 = 0.0;
  double t51 = 0.0;
  double t53 = 0.0;
  double t54 = 0.0;
  double t56 = 0.0;
  double t57 = 0.0;
  double t58 = 0.0;
  double t59 = 0.0;
  double t61 = 0.0;
  double t62 = 0.0;
  double t63 = 0.0;
  double t65 = 0.0;
  double t66 = 0.0;
  double t67 = 0.0;
  double t68 = 0.0;
  double t69 = 0.0;
  double t70 = 0.0;
  double t71 = 0.0;
  double t72 = 0.0;
  double t73 = 0.0;
  double t82 = 0.0;
  double t84 = 0.0;
  double t89 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; s = sigma[i % n_sigma];
    const double x1 = 1.0 + 1.0 / s; (void) x1;
    LG = R::lgammafn(x1);
    P0 = R::psigamma(x1, 0.0);
    P1 = R::psigamma(x1, 1.0);
    P2 = R::psigamma(x1, 2.0);
    P3 = R::psigamma(x1, 3.0);
    P4 = R::psigamma(x1, 4.0);
    t2 = LG*s;
    t3 = std::exp(t2);
    t6 = std::pow(s, 4);
    t8 = std::pow(s, 2);
    t10 = std::pow(s, 3);
    t15 = 6*LG;
    t17 = 4*LG;
    t18 = 12*P0;
    t19 = 6*s;
    t21 = std::pow(P0, 2);
    t22 = 3*t21;
    t23 = 6*P0;
    t26 = 2*s;
    t27 = std::pow(LG, 2);
    t28 = 1.0/s;
    t29 = P0*t28;
    t30 = 1.0/t8;
    t31 = P1*t30;
    t32 = 3*P1*t28;
    t34 = P0*t8;
    t37 = 3*t27;
    t38 = t37*t8;
    t39 = t21*t28;
    t41 = std::pow(LG, 3);
    t44 = LG*t22;
    t45 = 3*t31;
    t46 = 1.0/t6;
    t47 = 1.0/t10;
    t48 = 6*t39;
    t49 = std::pow(P0, 3);
    t51 = 4*t41;
    t53 = std::pow(LG, 4);
    t54 = std::pow(s, -6);
    t56 = 12*LG;
    t57 = std::pow(s, -5);
    t58 = 4*P2*t57;
    t59 = std::pow(P0, 4);
    t61 = P2*t46;
    t62 = LG*t18;
    t63 = std::pow(P1, 2);
    t65 = P1*t47;
    t66 = t21*t30;
    t67 = t30*t49;
    t68 = t18*t28;
    t69 = 6*t31;
    t70 = std::pow(s, -9);
    t71 = 60*P1;
    t72 = std::pow(s, -7);
    t73 = std::pow(s, -8);
    t82 = P3*t72;
    t84 = LG*t54;
    t89 = t46*t59;
  };
  const bool scalar = n_mean == 1 && n_sigma == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t0 = C*s;
    const double t1 = std::exp(t0);
    const double t4 = t1*t3;
    const double t5 = s*t4;
    const double t7 = t4*t6;
    const double t9 = t4*t8;
    const double t11 = t10*t4;
    const double t12 = 6*t4;
    const double t13 = 11*t9;
    const double t14 = 6*C;
    const double t16 = 4*C;
    const double t20 = 12*t0;
    const double t24 = t14*t8;
    const double t25 = std::pow(C, 2);
    const double t33 = 2*C;
    const double t35 = 3*t25;
    const double t36 = t35*t8;
    const double t40 = std::pow(C, 3);
    const double t42 = C*t22;
    const double t43 = t19*t25;
    const double t50 = 4*t40;
    const double t52 = std::pow(C, 4);
    const double t55 = 12*C;
    const double t60 = C*t18;
    const double t64 = LG*t55;
    const double t74 = 5*t4;
    const double t75 = C*t74;
    const double t76 = LG*t74;
    const double t77 = t4*t57;
    const double t78 = t71*t77;
    const double t79 = C*t78;
    const double t80 = C*t54;
    const double t81 = 40*P2*t4;
    const double t83 = LG*t78;
    const double t85 = 10*t4;
    const double t86 = t41*t85;
    const double t87 = t40*t85;
    const double t88 = 20*P0*P2*t4;
    const double t90 = 15*t4*t63;
    const double t91 = 30*t4;
    const double t92 = t65*t91;
    const double t93 = C*t27;
    const double t94 = 30*P1*t21*t77;
    const double t95 = LG*t25;
    const double t96 = t66*t91;
    o_mean_mean_mean_mean_mean[i] = s*(10*t11 + 24*t4 + 50*t5 + t7 + 35*t9 - 24)/std::pow(m, 5);
    o_mean_mean_mean_mean_sigma[i] = (-C*t13 - C*t7 - LG*t13 - LG*t7 + 11*P0*s*t1*t3 + P0*t1*t10*t3 + 6*P0*t1*t3*t8 + 6*P0*t1*t3 - t0*t12 - t11*t14 - t11*t15 - 4*t11 - t12*t2 - t12 - 22*t5 - 18*t9 + 6)/std::pow(m, 4);
    o_mean_mean_mean_sigma_sigma[i] = t4*(LG*t10*t33 + LG*t24 - 2*LG*t34 - P0*t16 - P0*t17 - P0*t19 + P1 + s*t21 + t0*t17 - t0*t23 + t10*t25 + t10*t27 + t15*t8 + t16 + t17 - t18 + t19 - t2*t23 + 12*t2 + t20 + t22 + t24 + t25*t26 + t26*t27 - 4*t29 + 2*t31 + t32 - t33*t34 + t36 + t38 + 2*t39 + 6)/std::pow(m, 3);
    o_mean_mean_sigma_sigma_sigma[i] = t4*(6*C*LG*P0*s + 6*C*LG*P0 + 6*C*P0*t28 + 12*C*P0 - C*t32 - C*t38 - C*t45 + 6*LG*P0*t28 + 12*LG*P0 - LG*t14 - LG*t20 - LG*t32 - LG*t36 - LG*t45 + 3*P0*P1*t30 + 3*P0*P1*t47 + 3*P0*s*t25 + 3*P0*s*t27 + 3*P0*t25 + 3*P0*t27 + 6*P0*t28 + P2*t46 + P2*t47 - s*t40 - s*t41 - t0*t37 - t14 - t15 - t19*t27 - t2*t35 - t22*t30 - t28*t42 - t28*t44 + t28*t49 + t30*t49 - t35 - t37 - t40*t8 - t41*t8 - t42 - t43 - t44 - t45 - t48)/std::pow(m, 2);
    o_mean_sigma_sigma_sigma_sigma[i] = t4*(-24*C*LG*t29 - P0*t50 - P0*t51 + P0*t58 + 6*P1*t21*t46 + P3*t54 + s*t52 + s*t53 + t0*t51 - t16*t61 - t16*t67 - t17*t61 - t17*t67 + t2*t50 + t25*t48 + t25*t56 - t25*t62 - t25*t68 + t25*t69 + t27*t43 + t27*t48 + t27*t55 - t27*t60 - t27*t68 + t27*t69 + t31*t64 + t39*t64 - 4*t47*t49 + t47*t59 + t50 + t51 + t55*t66 + t56*t66 + 3*t57*t63 + t58 - t60*t65 - t62*t65)/m;
    o_sigma_sigma_sigma_sigma_sigma[i] = -std::pow(C, 5)*t4 + 60*C*LG*P0*P1*t1*t3*t46 + 60*C*LG*P1*t1*t3*t46 + 20*C*LG*P2*t1*t3*t57 + 20*C*LG*t1*t3*t47*t49 + 20*C*P0*t1*t28*t3*t41 - C*t94 - std::pow(LG, 5)*t4 + 20*LG*P0*t1*t28*t3*t40 - LG*t94 + std::pow(P0, 5)*t1*t3*t57 + 30*P0*P1*t1*t25*t3*t46 + 30*P0*P1*t1*t27*t3*t46 + 60*P0*P1*t1*t3*t54 + 40*P0*P2*t1*t3*t72 + 5*P0*P3*t1*t3*t73 + 30*P0*t1*t25*t27*t28*t3 + 5*P0*t1*t28*t3*t52 + 5*P0*t1*t28*t3*t53 + 15*P0*t1*t3*t63*t72 - P0*t79 - P0*t83 + 10*P1*P2*t1*t3*t73 + 30*P1*t1*t21*t3*t54 + 30*P1*t1*t25*t3*t46 + 30*P1*t1*t27*t3*t46 + 10*P1*t1*t3*t49*t54 + 60*P1*t1*t3*t54 + 10*P2*t1*t21*t3*t72 + 10*P2*t1*t25*t3*t57 + 10*P2*t1*t27*t3*t57 + 60*P2*t1*t3*t72 - 60*P2*t72 + 15*P3*t1*t3*t73 - 15*P3*t73 + P4*t1*t3*t70 - P4*t70 + 10*t1*t25*t3*t47*t49 + 10*t1*t27*t3*t47*t49 + 30*t1*t3*t63*t72 - t25*t86 - t27*t87 - t52*t76 - t53*t75 - t54*t71 + 24*t57 - t65*t86 - t65*t87 - t66*t86 - t66*t87 - t75*t82 - t75*t89 - t76*t82 - t76*t89 - t79 - t80*t81 - t80*t88 - t80*t90 - t81*t84 - t83 - t84*t88 - t84*t90 - t92*t93 - t92*t95 - t93*t96 - t95*t96;
  }
  (void) LG;
  return List::create(Named("mean_mean_mean_mean_mean") = o_mean_mean_mean_mean_mean, Named("mean_mean_mean_mean_sigma") = o_mean_mean_mean_mean_sigma, Named("mean_mean_mean_sigma_sigma") = o_mean_mean_mean_sigma_sigma, Named("mean_mean_sigma_sigma_sigma") = o_mean_mean_sigma_sigma_sigma, Named("mean_sigma_sigma_sigma_sigma") = o_mean_sigma_sigma_sigma_sigma, Named("sigma_sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma_sigma);
}

// first derivatives of the expected Hessian
// [[Rcpp::export]]
List weibull3_dexpected1_cpp(NumericVector y, NumericVector mean, NumericVector sigma) {
  const int n = y.size();
  const int n_mean = mean.size(), n_sigma = sigma.size();
  NumericVector o_mean_mean_mean(n);
  NumericVector o_mean_mean_sigma(n);
  NumericVector o_sigma_sigma_mean(n);
  NumericVector o_sigma_sigma_sigma(n);
  NumericVector o_mean_sigma_mean(n);
  NumericVector o_mean_sigma_sigma(n);
  double m = 0.0, s = 0.0, LG = 0.0;
  double P0 = 0.0;
  double P1 = 0.0;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; s = sigma[i % n_sigma];
    const double x1 = 1.0 + 1.0 / s; (void) x1;
    LG = R::lgammafn(x1);
    P0 = R::psigamma(x1, 0.0);
    P1 = R::psigamma(x1, 1.0);
    t0 = std::pow(s, 2);
    t1 = std::pow(m, -2);
    t2 = P1/s;
  };
  const bool scalar = n_mean == 1 && n_sigma == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    o_mean_mean_mean[i] = 2.0*t0/std::pow(m, 3);
    o_mean_mean_sigma[i] = -2.0*s*t1;
    o_sigma_sigma_mean[i] = 0;
    o_sigma_sigma_sigma[i] = (2.0*std::pow(P0, 2) + 2.0*P0*t2 - 1.6911373403938686*P0 - 0.84556867019693428*t2 + 3.6473613217057588)/std::pow(s, 3);
    o_mean_sigma_mean[i] = t1*(1.0*P0 - 0.42278433509846714);
    o_mean_sigma_sigma[i] = 1.0*P1/(m*t0);
  }
  (void) LG;
  return List::create(Named("mean_mean_mean") = o_mean_mean_mean, Named("mean_mean_sigma") = o_mean_mean_sigma, Named("sigma_sigma_mean") = o_sigma_sigma_mean, Named("sigma_sigma_sigma") = o_sigma_sigma_sigma, Named("mean_sigma_mean") = o_mean_sigma_mean, Named("mean_sigma_sigma") = o_mean_sigma_sigma);
}

// second derivatives of the expected Hessian
// [[Rcpp::export]]
List weibull3_dexpected2_cpp(NumericVector y, NumericVector mean, NumericVector sigma) {
  const int n = y.size();
  const int n_mean = mean.size(), n_sigma = sigma.size();
  NumericVector o_mean_mean_mean_mean(n);
  NumericVector o_mean_mean_sigma_sigma(n);
  NumericVector o_mean_mean_mean_sigma(n);
  NumericVector o_sigma_sigma_mean_mean(n);
  NumericVector o_sigma_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_mean_sigma(n);
  NumericVector o_mean_sigma_mean_mean(n);
  NumericVector o_mean_sigma_sigma_sigma(n);
  NumericVector o_mean_sigma_mean_sigma(n);
  double m = 0.0, s = 0.0, LG = 0.0;
  double P0 = 0.0;
  double P1 = 0.0;
  double P2 = 0.0;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double t3 = 0.0;
  double t4 = 0.0;
  double t5 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; s = sigma[i % n_sigma];
    const double x1 = 1.0 + 1.0 / s; (void) x1;
    LG = R::lgammafn(x1);
    P0 = R::psigamma(x1, 0.0);
    P1 = R::psigamma(x1, 1.0);
    P2 = R::psigamma(x1, 2.0);
    t0 = std::pow(s, 2);
    t1 = std::pow(m, -2);
    t2 = std::pow(m, -3);
    t3 = 1.0/s;
    t4 = 1.0/t0;
    t5 = 2.0*P0;
  };
  const bool scalar = n_mean == 1 && n_sigma == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    o_mean_mean_mean_mean[i] = -6.0*t0/std::pow(m, 4);
    o_mean_mean_sigma_sigma[i] = -2.0*t1;
    o_mean_mean_mean_sigma[i] = 4.0*s*t2;
    o_sigma_sigma_mean_mean[i] = 0;
    o_sigma_sigma_sigma_sigma[i] = (-6.0*std::pow(P0, 2) - 12.0*P0*P1*t3 + 5.0734120211816057*P0 - 2.0*std::pow(P1, 2)*t4 + 5.0734120211816057*P1*t3 - P2*t4*t5 + 0.84556867019693428*P2*t4 - 10.942083965117276)/std::pow(s, 4);
    o_sigma_sigma_mean_sigma[i] = 0;
    o_mean_sigma_mean_mean[i] = t2*(0.84556867019693428 - t5);
    o_mean_sigma_sigma_sigma[i] = -(2.0*P1 + 1.0*P2*t3)/(m*std::pow(s, 3));
    o_mean_sigma_mean_sigma[i] = -1.0*P1*t1*t4;
  }
  (void) LG;
  return List::create(Named("mean_mean_mean_mean") = o_mean_mean_mean_mean, Named("mean_mean_sigma_sigma") = o_mean_mean_sigma_sigma, Named("mean_mean_mean_sigma") = o_mean_mean_mean_sigma, Named("sigma_sigma_mean_mean") = o_sigma_sigma_mean_mean, Named("sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma, Named("sigma_sigma_mean_sigma") = o_sigma_sigma_mean_sigma, Named("mean_sigma_mean_mean") = o_mean_sigma_mean_mean, Named("mean_sigma_sigma_sigma") = o_mean_sigma_sigma_sigma, Named("mean_sigma_mean_sigma") = o_mean_sigma_mean_sigma);
}

// order 1 of the distribution function in (mean, sigma)
// [[Rcpp::export]]
List weibull3_dcdf1_cpp(NumericVector y, NumericVector mean, NumericVector sigma) {
  const int n = y.size();
  const int n_mean = mean.size(), n_sigma = sigma.size();
  NumericVector o_mean(n);
  NumericVector o_sigma(n);
  double m = 0.0, s = 0.0, LG = 0.0;
  double P0 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; s = sigma[i % n_sigma];
    const double x1 = 1.0 + 1.0 / s; (void) x1;
    LG = R::lgammafn(x1);
    P0 = R::psigamma(x1, 0.0);
  };
  const bool scalar = n_mean == 1 && n_sigma == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t0 = std::exp(C*s)*std::exp(LG*s);
    const double t1 = t0*std::exp(-t0);
    o_mean[i] = -s*t1/m;
    o_sigma[i] = t1*(C + LG - P0/s);
  }
  (void) LG;
  return List::create(Named("mean") = o_mean, Named("sigma") = o_sigma);
}

// order 2 of the distribution function in (mean, sigma)
// [[Rcpp::export]]
List weibull3_dcdf2_cpp(NumericVector y, NumericVector mean, NumericVector sigma) {
  const int n = y.size();
  const int n_mean = mean.size(), n_sigma = sigma.size();
  NumericVector o_mean_mean(n);
  NumericVector o_sigma_sigma(n);
  NumericVector o_mean_sigma(n);
  double m = 0.0, s = 0.0, LG = 0.0;
  double P0 = 0.0;
  double P1 = 0.0;
  double t2 = 0.0;
  double t3 = 0.0;
  double t8 = 0.0;
  double t11 = 0.0;
  double t12 = 0.0;
  double t13 = 0.0;
  double t14 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; s = sigma[i % n_sigma];
    const double x1 = 1.0 + 1.0 / s; (void) x1;
    LG = R::lgammafn(x1);
    P0 = R::psigamma(x1, 0.0);
    P1 = R::psigamma(x1, 1.0);
    t2 = LG*s;
    t3 = std::exp(t2);
    t8 = std::pow(LG, 2);
    t11 = std::pow(P0, 2)/std::pow(s, 2);
    t12 = 1.0/s;
    t13 = P0*t12;
    t14 = 2*LG;
  };
  const bool scalar = n_mean == 1 && n_sigma == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t0 = C*s;
    const double t1 = std::exp(t0);
    const double t4 = t1*t3;
    const double t5 = std::exp(-t4);
    const double t6 = s*t4;
    const double t7 = std::pow(C, 2);
    const double t9 = 2*C;
    const double t10 = LG*t9;
    const double t15 = P0*t4;
    const double t16 = t12*t15;
    const double t17 = t4*t5;
    o_mean_mean[i] = t5*t6*(s - t6 + 1)/std::pow(m, 2);
    o_sigma_sigma[i] = t17*(P1/std::pow(s, 3) - t10*t4 + t10 - t11*t4 + t11 - t13*t14 - t13*t9 + t14*t16 + t16*t9 - t4*t7 - t4*t8 + t7 + t8);
    o_mean_sigma[i] = t17*(C*s*t1*t3 + LG*s*t1*t3 + P0 - t0 - t15 - t2 - 1)/m;
  }
  (void) LG;
  return List::create(Named("mean_mean") = o_mean_mean, Named("sigma_sigma") = o_sigma_sigma, Named("mean_sigma") = o_mean_sigma);
}

// order 3 of the distribution function in (mean, sigma)
// [[Rcpp::export]]
List weibull3_dcdf3_cpp(NumericVector y, NumericVector mean, NumericVector sigma) {
  const int n = y.size();
  const int n_mean = mean.size(), n_sigma = sigma.size();
  NumericVector o_mean_mean_mean(n);
  NumericVector o_mean_mean_sigma(n);
  NumericVector o_mean_sigma_sigma(n);
  NumericVector o_sigma_sigma_sigma(n);
  double m = 0.0, s = 0.0, LG = 0.0;
  double P0 = 0.0;
  double P1 = 0.0;
  double P2 = 0.0;
  double t0 = 0.0;
  double t1 = 0.0;
  double t4 = 0.0;
  double t5 = 0.0;
  double t7 = 0.0;
  double t11 = 0.0;
  double t12 = 0.0;
  double t14 = 0.0;
  double t17 = 0.0;
  double t20 = 0.0;
  double t21 = 0.0;
  double t22 = 0.0;
  double t24 = 0.0;
  double t25 = 0.0;
  double t26 = 0.0;
  double t30 = 0.0;
  double t31 = 0.0;
  double t35 = 0.0;
  double t36 = 0.0;
  double t37 = 0.0;
  double t38 = 0.0;
  double t40 = 0.0;
  double t41 = 0.0;
  double t42 = 0.0;
  double t44 = 0.0;
  double t45 = 0.0;
  double t48 = 0.0;
  double t49 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; s = sigma[i % n_sigma];
    const double x1 = 1.0 + 1.0 / s; (void) x1;
    LG = R::lgammafn(x1);
    P0 = R::psigamma(x1, 0.0);
    P1 = R::psigamma(x1, 1.0);
    P2 = R::psigamma(x1, 2.0);
    t0 = 3*s;
    t1 = std::pow(s, 2);
    t4 = LG*s;
    t5 = std::exp(t4);
    t7 = std::exp(2*t4);
    t11 = 2*s;
    t12 = P0*s;
    t14 = LG*t1;
    t17 = 2*LG;
    t20 = std::pow(LG, 2);
    t21 = s*t20;
    t22 = 1.0/t1;
    t24 = 1.0/s;
    t25 = std::pow(P0, 2);
    t26 = t24*t25;
    t30 = std::pow(LG, 3);
    t31 = 3*t20;
    t35 = 3*P1;
    t36 = t35/std::pow(s, 4);
    t37 = std::pow(s, -3);
    t38 = t35*t37;
    t40 = LG*t38;
    t41 = std::pow(P0, 3)*t37;
    t42 = P0*t24;
    t44 = t22*t25;
    t45 = 3*t44;
    t48 = LG*t45;
    t49 = t31*t42;
  };
  const bool scalar = n_mean == 1 && n_sigma == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t2 = C*s;
    const double t3 = std::exp(t2);
    const double t6 = std::exp(2*t2);
    const double t8 = t6*t7;
    const double t9 = t3*t5;
    const double t10 = t9*std::exp(-t9);
    const double t13 = C*t1;
    const double t15 = P0*t9;
    const double t16 = 3*t9;
    const double t18 = std::pow(C, 2);
    const double t19 = s*t18;
    const double t23 = t17*t2;
    const double t27 = 6*t15;
    const double t28 = t15*t24;
    const double t29 = std::pow(C, 3);
    const double t32 = C*t31;
    const double t33 = 3*t18;
    const double t34 = LG*t33;
    const double t39 = C*t38;
    const double t43 = 6*C*LG*t42;
    const double t46 = C*t45;
    const double t47 = t33*t42;
    const double t50 = 9*t9;
    const double t51 = C*t50;
    const double t52 = LG*t50;
    const double t53 = 9*t28;
    o_mean_mean_mean[i] = s*t10*(3*s*t3*t5 - t0 + 3*t1*t3*t5 - t1*t8 - t1 - 2)/std::pow(m, 3);
    o_mean_mean_sigma[i] = t10*(-P0 + t0*t15 - t11*t9 + t11 - t12*t8 - t12 - t13*t16 + t13*t8 + t13 - t14*t16 + t14*t8 + t14 + t15 - t2*t9 + t2 - t4*t9 + t4 + 1)/std::pow(m, 2);
    o_mean_sigma_sigma[i] = t10*(6*C*LG*s*t3*t5 + 2*C*P0*t6*t7 + 2*C*P0 - C*t27 + 2*C*t3*t5 - 2*C + 2*LG*P0*t6*t7 + 2*LG*P0 - LG*t27 + 2*LG*t3*t5 + 2*P0*t24 + P1*t22*t3*t5 - P1*t22 + 3*s*t18*t3*t5 + 3*s*t20*t3*t5 - t17 - t19*t8 - t19 - t21*t8 - t21 - t23*t8 - t23 + 3*t24*t25*t3*t5 - t26*t8 - t26 - 2*t28)/m;
    o_sigma_sigma_sigma[i] = t10*(18*C*LG*t28 - P0*t36 - P2/std::pow(s, 5) + t15*t36 - t16*t29 - t16*t30 + t16*t41 - t18*t52 + t18*t53 - t20*t51 + t20*t53 + t29*t8 + t29 + t30*t8 + t30 + t32*t8 + t32 + t34*t8 + t34 - t36 - t39*t9 + t39 - t40*t9 + t40 - t41*t8 - t41 - t43*t8 - t43 - t44*t51 - t44*t52 + t46*t8 + t46 - t47*t8 - t47 + t48*t8 + t48 - t49*t8 - t49);
  }
  (void) LG;
  return List::create(Named("mean_mean_mean") = o_mean_mean_mean, Named("mean_mean_sigma") = o_mean_mean_sigma, Named("mean_sigma_sigma") = o_mean_sigma_sigma, Named("sigma_sigma_sigma") = o_sigma_sigma_sigma);
}

// order 4 of the distribution function in (mean, sigma)
// [[Rcpp::export]]
List weibull3_dcdf4_cpp(NumericVector y, NumericVector mean, NumericVector sigma) {
  const int n = y.size();
  const int n_mean = mean.size(), n_sigma = sigma.size();
  NumericVector o_mean_mean_mean_mean(n);
  NumericVector o_mean_mean_mean_sigma(n);
  NumericVector o_mean_mean_sigma_sigma(n);
  NumericVector o_mean_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_sigma_sigma(n);
  double m = 0.0, s = 0.0, LG = 0.0;
  double P0 = 0.0;
  double P1 = 0.0;
  double P2 = 0.0;
  double P3 = 0.0;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  double t3 = 0.0;
  double t6 = 0.0;
  double t7 = 0.0;
  double t13 = 0.0;
  double t14 = 0.0;
  double t19 = 0.0;
  double t22 = 0.0;
  double t23 = 0.0;
  double t25 = 0.0;
  double t27 = 0.0;
  double t28 = 0.0;
  double t31 = 0.0;
  double t33 = 0.0;
  double t34 = 0.0;
  double t35 = 0.0;
  double t38 = 0.0;
  double t39 = 0.0;
  double t42 = 0.0;
  double t43 = 0.0;
  double t44 = 0.0;
  double t45 = 0.0;
  double t46 = 0.0;
  double t47 = 0.0;
  double t50 = 0.0;
  double t51 = 0.0;
  double t53 = 0.0;
  double t54 = 0.0;
  double t55 = 0.0;
  double t59 = 0.0;
  double t70 = 0.0;
  double t72 = 0.0;
  double t73 = 0.0;
  double t76 = 0.0;
  double t77 = 0.0;
  double t78 = 0.0;
  double t79 = 0.0;
  double t80 = 0.0;
  double t81 = 0.0;
  double t85 = 0.0;
  double t86 = 0.0;
  double t88 = 0.0;
  double t92 = 0.0;
  double t95 = 0.0;
  double t97 = 0.0;
  double t98 = 0.0;
  double t99 = 0.0;
  double t100 = 0.0;
  double t101 = 0.0;
  double t102 = 0.0;
  double t103 = 0.0;
  double t105 = 0.0;
  double t107 = 0.0;
  double t108 = 0.0;
  double t109 = 0.0;
  double t112 = 0.0;
  double t115 = 0.0;
  double t116 = 0.0;
  double t119 = 0.0;
  double t121 = 0.0;
  double t122 = 0.0;
  double t123 = 0.0;
  double t124 = 0.0;
  double t125 = 0.0;
  double t127 = 0.0;
  double t128 = 0.0;
  double t143 = 0.0;
  double t145 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; s = sigma[i % n_sigma];
    const double x1 = 1.0 + 1.0 / s; (void) x1;
    LG = R::lgammafn(x1);
    P0 = R::psigamma(x1, 0.0);
    P1 = R::psigamma(x1, 1.0);
    P2 = R::psigamma(x1, 2.0);
    P3 = R::psigamma(x1, 3.0);
    t0 = 11*s;
    t1 = std::pow(s, 3);
    t2 = std::pow(s, 2);
    t3 = 6*t2;
    t6 = LG*s;
    t7 = std::exp(t6);
    t13 = 2*t6;
    t14 = std::exp(t13);
    t19 = std::exp(3*t6);
    t22 = 6*s;
    t23 = 3*t2;
    t25 = LG*t1;
    t27 = LG*t23;
    t28 = 2*P0;
    t31 = P0*t2;
    t33 = 2*LG;
    t34 = 4*P0;
    t35 = std::pow(P0, 2);
    t38 = LG*t28;
    t39 = 4*t6;
    t42 = std::pow(LG, 2);
    t43 = s*t42;
    t44 = 1.0/t2;
    t45 = P1*t44;
    t46 = 1.0/s;
    t47 = P1*t46;
    t50 = t28*t6;
    t51 = t28*t46;
    t53 = t2*t42;
    t54 = t35*t46;
    t55 = LG*t2;
    t59 = 6*LG;
    t70 = 3*t42;
    t72 = std::pow(LG, 3);
    t73 = std::pow(s, -4);
    t76 = 3*LG;
    t77 = t45*t76;
    t78 = 1.0/t1;
    t79 = t35*t44;
    t80 = 3*t79;
    t81 = std::pow(P0, 3);
    t85 = 18*LG;
    t86 = P1*t78;
    t88 = t44*t81;
    t92 = std::pow(LG, 4);
    t95 = 4*LG;
    t97 = std::pow(s, -5);
    t98 = 12*P1;
    t99 = t97*t98;
    t100 = std::pow(s, -6);
    t101 = P2*t100;
    t102 = std::pow(P0, 4)*t73;
    t103 = t73*t98;
    t105 = P2*t97;
    t107 = LG*t103;
    t108 = t105*t95;
    t109 = t101*t34;
    t112 = 3*std::pow(P1, 2)*t100;
    t115 = P0*t107;
    t116 = t78*t81;
    t119 = t34*t46;
    t121 = t116*t95;
    t122 = 6*t42;
    t123 = t119*t72;
    t124 = P1*t97;
    t125 = t124*t35;
    t127 = P0*t46;
    t128 = 12*t127;
    t143 = 18*t86;
    t145 = LG*t116;
  };
  const bool scalar = n_mean == 1 && n_sigma == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t4 = C*s;
    const double t5 = std::exp(t4);
    const double t8 = t5*t7;
    const double t9 = t2*t8;
    const double t10 = 7*t8;
    const double t11 = 2*t4;
    const double t12 = std::exp(t11);
    const double t15 = t12*t14;
    const double t16 = t15*t3;
    const double t17 = 6*t15;
    const double t18 = std::exp(3*t4);
    const double t20 = t18*t19;
    const double t21 = t8*std::exp(-t8);
    const double t24 = C*t1;
    const double t26 = C*t23;
    const double t29 = P0*t8;
    const double t30 = 9*t29;
    const double t32 = 2*C;
    const double t36 = C*t28;
    const double t37 = 4*t4;
    const double t40 = std::pow(C, 2);
    const double t41 = s*t40;
    const double t48 = t33*t4;
    const double t49 = t28*t4;
    const double t52 = t2*t40;
    const double t56 = t32*t55;
    const double t57 = 6*C;
    const double t58 = 12*t8;
    const double t60 = 14*t29;
    const double t61 = 3*t8;
    const double t62 = t17*t35;
    const double t63 = t20*t35;
    const double t64 = C*LG;
    const double t65 = 12*t15;
    const double t66 = P0*t65;
    const double t67 = t15*t54;
    const double t68 = LG*t57;
    const double t69 = 3*t40;
    const double t71 = std::pow(C, 3);
    const double t74 = 3*C;
    const double t75 = t45*t74;
    const double t82 = 21*t29;
    const double t83 = 18*C;
    const double t84 = t29*t46;
    const double t87 = t15*t22;
    const double t89 = t15*t42;
    const double t90 = t15*t40;
    const double t91 = std::pow(C, 4);
    const double t93 = 4*C;
    const double t94 = t72*t93;
    const double t96 = t71*t95;
    const double t104 = C*t103;
    const double t106 = t105*t93;
    const double t110 = 6*t40;
    const double t111 = t110*t42;
    const double t113 = t64*t78*t98;
    const double t114 = P0*t104;
    const double t117 = t116*t93;
    const double t118 = t110*t86;
    const double t120 = t119*t71;
    const double t126 = 12*t64;
    const double t129 = C*t42;
    const double t130 = t128*t129;
    const double t131 = LG*t40;
    const double t132 = t128*t131;
    const double t133 = C*t8;
    const double t134 = 28*t133;
    const double t135 = LG*t71;
    const double t136 = 28*t8;
    const double t137 = 24*t15;
    const double t138 = C*t137;
    const double t139 = t42*t8;
    const double t140 = 42*t40;
    const double t141 = LG*t133;
    const double t142 = 36*P1*t29*t73;
    const double t144 = 28*t84;
    const double t146 = 36*t90;
    const double t147 = 84*t84;
    const double t148 = t127*t137;
    const double t149 = t44*t63;
    o_mean_mean_mean_mean[i] = s*t21*(-t0*t8 + t0 - t1*t10 + t1*t17 - t1*t20 + t1 + t16 + t3 - 18*t9 + 6)/std::pow(m, 4);
    o_mean_mean_mean_sigma[i] = t21*(2*C*s*t5*t7 + C*t1*t18*t19 + 7*C*t1*t5*t7 + 9*C*t2*t5*t7 + 2*LG*s*t5*t7 + LG*t1*t18*t19 + 7*LG*t1*t5*t7 + 9*LG*t2*t5*t7 + 3*P0*s*t12*t14 + 3*P0*s + 6*P0*t12*t14*t2 + P0*t2 + 2*P0 - s*t30 + 6*s*t5*t7 - t10*t31 - t11 - t13 - t15*t23 - t15*t26 - t15*t27 - t17*t24 - t17*t25 + 9*t2*t5*t7 - t20*t31 - t22 - t23 - t24 - t25 - t26 - t27 - t28*t8 - 2)/std::pow(m, 3);
    o_mean_mean_sigma_sigma[i] = t21*(C*t55*t65 - t10*t35 - t10*t52 - t10*t53 - t15*t34 - t15*t36 + t15*t37 - t15*t38 + t15*t39 + t15*t41 + t15*t43 + t15*t47 + t15*t48 + t16*t40 + t16*t42 + t20*t49 + t20*t50 - t20*t52 - t20*t53 - t20*t56 + t29*t57 + t29*t59 + 12*t29 - t32*t8 + t32 - t33*t8 + t33 - t34 + t35 - t36 + t37 - t38 + t39 - t4*t58 - t4*t59*t8 + t4*t60 - t4*t66 - t41*t61 + t41 - t43*t61 + t43 - t45*t8 + t45 - t47*t61 + t47 + t48 - t49 - t50 + t51*t8 - t51 + t52 + t53 - t54*t61 + t54 + t56 - t58*t6 + t6*t60 - t6*t66 + t62 - t63 - 14*t64*t9 + t67 - 2*t8 + 2)/std::pow(m, 2);
    o_mean_sigma_sigma_sigma[i] = t21*(36*C*LG*P0*t12*t14 + 6*C*LG*P0 + 18*C*LG*t5*t7 + 6*C*P0*t12*t14*t46 + 6*C*P0*t46 + 9*C*P1*t44*t5*t7 + 3*C*s*t18*t19*t42 + 21*C*s*t42*t5*t7 + 3*C*t18*t19*t35*t46 + 21*C*t35*t46*t5*t7 + 6*LG*P0*t12*t14*t46 + 6*LG*P0*t46 + 9*LG*P1*t44*t5*t7 + 3*LG*s*t18*t19*t40 + 21*LG*s*t40*t5*t7 + 3*LG*t18*t19*t35*t46 + 21*LG*t35*t46*t5*t7 + 3*P0*P1*t12*t14*t78 + 3*P0*P1*t78 + 18*P0*t12*t14*t40 + 18*P0*t12*t14*t42 - P0*t20*t68 - P0*t20*t69 - P0*t20*t70 + 3*P0*t40 + 3*P0*t42 - P2*t73*t8 + P2*t73 + s*t18*t19*t71 + s*t18*t19*t72 + 7*s*t5*t7*t71 + 7*s*t5*t7*t72 - s*t71 - s*t72 - t10*t88 + 6*t12*t14*t44*t81 - t15*t68 - t15*t69 - t15*t70 - t15*t75 - t15*t77 - t15*t80 - t20*t88 - 42*t29*t64 - t30*t86 + 9*t35*t44*t5*t7 - t4*t70 - 18*t4*t89 + 9*t40*t5*t7 - t40*t82 + 9*t42*t5*t7 - t42*t82 + t44*t81 - t54*t74 - t54*t76 - t6*t69 - 18*t6*t90 - t67*t83 - t67*t85 - t68 - t69 - t70 - t71*t87 - t72*t87 - t75 - t77 - t80 - t83*t84 - t84*t85)/m;
    o_sigma_sigma_sigma_sigma[i] = t21*(-72*C*t127*t89 + C*t142 - 72*LG*t127*t90 + LG*t142 + P0*t99 + P3/std::pow(s, 7) - t10*t102 - t10*t91 - t10*t92 + 8*t101 + t102*t17 - t102*t20 + t102 + t104*t8 - t104 + t106*t8 - t106 + t107*t8 - t107 + t108*t8 - t108 - t109*t8 + t109 - t110*t149 + t110*t79 - t111*t20 + t111 - t112*t8 + t112 + t113*t15 + t113 - t114*t15 - t114 - t115*t15 - t115 + t116*t134 - t116*t138 + t117*t20 - t117 + t118*t15 + t118 + t120*t20 - t120 + t121*t20 - t121 - t122*t149 + t122*t79 + t122*t86 + t123*t20 - t123 + t124*t62 - 18*t125*t8 + 6*t125 - t126*t149 + t126*t79 + t129*t147 + t130*t20 - t130 + t131*t147 + t132*t20 - t132 - t134*t72 - t135*t136 + t135*t137 + t136*t145 - t137*t145 + t138*t72 - t139*t140 - t139*t143 - 42*t139*t79 - t140*t79*t8 - 84*t141*t79 - 36*t141*t86 - t143*t40*t8 + t144*t71 + t144*t72 + t146*t42 + t146*t79 - t148*t71 - t148*t72 + 72*t15*t64*t79 + t17*t42*t86 + t17*t91 + t17*t92 - t20*t91 - t20*t92 - t20*t94 - t20*t96 - t29*t99 + 36*t79*t89 + t91 + t92 + t94 + t96 + t99);
  }
  (void) LG;
  return List::create(Named("mean_mean_mean_mean") = o_mean_mean_mean_mean, Named("mean_mean_mean_sigma") = o_mean_mean_mean_sigma, Named("mean_mean_sigma_sigma") = o_mean_mean_sigma_sigma, Named("mean_sigma_sigma_sigma") = o_mean_sigma_sigma_sigma, Named("sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma);
}

// 1 in y and 1 in (mean, sigma)
// [[Rcpp::export]]
List weibull3_cross_y_cpp(NumericVector y, NumericVector mean, NumericVector sigma) {
  const int n = y.size();
  const int n_mean = mean.size(), n_sigma = sigma.size();
  NumericVector o_mean(n);
  NumericVector o_sigma(n);
  double m = 0.0, s = 0.0, LG = 0.0;
  double P0 = 0.0;
  double t3 = 0.0;
  double t4 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; s = sigma[i % n_sigma];
    const double x1 = 1.0 + 1.0 / s; (void) x1;
    LG = R::lgammafn(x1);
    P0 = R::psigamma(x1, 0.0);
    t3 = LG*s;
    t4 = std::exp(t3);
  };
  const bool scalar = n_mean == 1 && n_sigma == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t0 = 1.0/yy;
    const double t1 = C*s;
    const double t2 = std::exp(t1);
    const double t5 = t2*t4;
    o_mean[i] = std::pow(s, 2)*t0*t5/m;
    o_sigma[i] = t0*(P0*t2*t4 - t1*t5 - t3*t5 - t5 + 1);
  }
  (void) LG;
  return List::create(Named("mean") = o_mean, Named("sigma") = o_sigma);
}

// 2 in y and 1 in (mean, sigma)
// [[Rcpp::export]]
List weibull3_cross2_y_cpp(NumericVector y, NumericVector mean, NumericVector sigma) {
  const int n = y.size();
  const int n_mean = mean.size(), n_sigma = sigma.size();
  NumericVector o_mean(n);
  NumericVector o_sigma(n);
  double m = 0.0, s = 0.0, LG = 0.0;
  double P0 = 0.0;
  double t2 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; s = sigma[i % n_sigma];
    const double x1 = 1.0 + 1.0 / s; (void) x1;
    LG = R::lgammafn(x1);
    P0 = R::psigamma(x1, 0.0);
    t2 = std::exp(LG*s);
  };
  const bool scalar = n_mean == 1 && n_sigma == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t0 = std::pow(yy, -2);
    const double t1 = std::exp(C*s);
    const double t3 = t1*t2;
    const double t4 = std::pow(s, 2)*t3;
    o_mean[i] = t0*t4*(s - 1)/m;
    o_sigma[i] = t0*(C*s*t1*t2 - C*t4 + LG*s*t1*t2 - LG*t4 + P0*s*t1*t2 - P0*t3 - 2*s*t3 + t1*t2 - 1);
  }
  (void) LG;
  return List::create(Named("mean") = o_mean, Named("sigma") = o_sigma);
}

// 1 in y and 2 in (mean, sigma)
// [[Rcpp::export]]
List weibull3_grad_y_hess_cpp(NumericVector y, NumericVector mean, NumericVector sigma) {
  const int n = y.size();
  const int n_mean = mean.size(), n_sigma = sigma.size();
  NumericVector o_mean_mean(n);
  NumericVector o_sigma_sigma(n);
  NumericVector o_mean_sigma(n);
  double m = 0.0, s = 0.0, LG = 0.0;
  double P0 = 0.0;
  double P1 = 0.0;
  double t0 = 0.0;
  double t2 = 0.0;
  double t4 = 0.0;
  double t5 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; s = sigma[i % n_sigma];
    const double x1 = 1.0 + 1.0 / s; (void) x1;
    LG = R::lgammafn(x1);
    P0 = R::psigamma(x1, 0.0);
    P1 = R::psigamma(x1, 1.0);
    t0 = std::pow(s, 2);
    t2 = LG*s;
    t4 = 2*LG;
    t5 = 1.0/s;
  };
  const bool scalar = n_mean == 1 && n_sigma == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t1 = C*s;
    const double t3 = std::exp(t1)*std::exp(t2)/yy;
    o_mean_mean[i] = -t0*t3*(s + 1)/std::pow(m, 2);
    o_sigma_sigma[i] = t3*(-std::pow(C, 2)*s + 2*C*P0 - 2*C - std::pow(LG, 2)*s + 2*LG*P0 - std::pow(P0, 2)*t5 + 2*P0*t5 - P1/t0 - t1*t4 - t4);
    o_mean_sigma[i] = s*t3*(-P0 + t1 + t2 + 2)/m;
  }
  (void) LG;
  return List::create(Named("mean_mean") = o_mean_mean, Named("sigma_sigma") = o_sigma_sigma, Named("mean_sigma") = o_mean_sigma);
}

// 2 in y and 2 in (mean, sigma)
// [[Rcpp::export]]
List weibull3_hess_y_hess_cpp(NumericVector y, NumericVector mean, NumericVector sigma) {
  const int n = y.size();
  const int n_mean = mean.size(), n_sigma = sigma.size();
  NumericVector o_mean_mean(n);
  NumericVector o_sigma_sigma(n);
  NumericVector o_mean_sigma(n);
  double m = 0.0, s = 0.0, LG = 0.0;
  double P0 = 0.0;
  double P1 = 0.0;
  double t0 = 0.0;
  double t2 = 0.0;
  double t4 = 0.0;
  double t6 = 0.0;
  double t8 = 0.0;
  double t9 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; s = sigma[i % n_sigma];
    const double x1 = 1.0 + 1.0 / s; (void) x1;
    LG = R::lgammafn(x1);
    P0 = R::psigamma(x1, 0.0);
    P1 = R::psigamma(x1, 1.0);
    t0 = std::pow(s, 2);
    t2 = LG*s;
    t4 = std::pow(P0, 2);
    t6 = 1.0/s;
    t8 = std::pow(LG, 2);
    t9 = LG*t0;
  };
  const bool scalar = n_mean == 1 && n_sigma == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t1 = C*s;
    const double t3 = std::exp(t1)*std::exp(t2)/std::pow(yy, 2);
    const double t5 = 2*C;
    const double t7 = std::pow(C, 2);
    o_mean_mean[i] = t0*t3*(1 - t0)/std::pow(m, 2);
    o_sigma_sigma[i] = t3*(2*C*LG*s + 2*C*P0*s + 2*C + 2*LG*P0*s - 2*LG*P0 + 2*LG - P0*t5 - 2*P0*t6 + 4*P0 - P1*t6 + P1/t0 + s*t7 + s*t8 - t0*t7 - t0*t8 - 4*t1 - 4*t2 + t4*t6 - t4 - t5*t9 - 2);
    o_mean_sigma[i] = s*t3*(C*t0 - P0*s + P0 + 3*s - t1 - t2 + t9 - 2)/m;
  }
  (void) LG;
  return List::create(Named("mean_mean") = o_mean_mean, Named("sigma_sigma") = o_sigma_sigma, Named("mean_sigma") = o_mean_sigma);
}

// order 1 of log f in y
// [[Rcpp::export]]
List weibull3_dy1_cpp(NumericVector y, NumericVector mean, NumericVector sigma) {
  const int n = y.size();
  const int n_mean = mean.size(), n_sigma = sigma.size();
  NumericVector o_y(n);
  double m = 0.0, s = 0.0, LG = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; s = sigma[i % n_sigma];
    const double x1 = 1.0 + 1.0 / s; (void) x1;
    LG = R::lgammafn(x1);
  };
  const bool scalar = n_mean == 1 && n_sigma == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    o_y[i] = (-s*std::exp(C*s)*std::exp(LG*s) + s - 1)/yy;
  }
  (void) LG;
  return List::create(Named("y") = o_y);
}

// order 2 of log f in y
// [[Rcpp::export]]
List weibull3_dy2_cpp(NumericVector y, NumericVector mean, NumericVector sigma) {
  const int n = y.size();
  const int n_mean = mean.size(), n_sigma = sigma.size();
  NumericVector o_y(n);
  double m = 0.0, s = 0.0, LG = 0.0;
  double t1 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; s = sigma[i % n_sigma];
    const double x1 = 1.0 + 1.0 / s; (void) x1;
    LG = R::lgammafn(x1);
    t1 = std::exp(LG*s);
  };
  const bool scalar = n_mean == 1 && n_sigma == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t0 = std::exp(C*s);
    o_y[i] = (-std::pow(s, 2)*t0*t1 + s*t0*t1 - s + 1)/std::pow(yy, 2);
  }
  (void) LG;
  return List::create(Named("y") = o_y);
}

// order 3 of log f in y
// [[Rcpp::export]]
List weibull3_dy3_cpp(NumericVector y, NumericVector mean, NumericVector sigma) {
  const int n = y.size();
  const int n_mean = mean.size(), n_sigma = sigma.size();
  NumericVector o_y(n);
  double m = 0.0, s = 0.0, LG = 0.0;
  double t1 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; s = sigma[i % n_sigma];
    const double x1 = 1.0 + 1.0 / s; (void) x1;
    LG = R::lgammafn(x1);
    t1 = std::exp(LG*s);
  };
  const bool scalar = n_mean == 1 && n_sigma == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t0 = std::exp(C*s);
    const double t2 = t0*t1;
    o_y[i] = (-std::pow(s, 3)*t2 + 3*std::pow(s, 2)*t0*t1 - 2*s*t2 + 2*s - 2)/std::pow(yy, 3);
  }
  (void) LG;
  return List::create(Named("y") = o_y);
}

// order 4 of log f in y
// [[Rcpp::export]]
List weibull3_dy4_cpp(NumericVector y, NumericVector mean, NumericVector sigma) {
  const int n = y.size();
  const int n_mean = mean.size(), n_sigma = sigma.size();
  NumericVector o_y(n);
  double m = 0.0, s = 0.0, LG = 0.0;
  double t1 = 0.0;
  auto set_params = [&](int i) {
    m = mean[i % n_mean]; s = sigma[i % n_sigma];
    const double x1 = 1.0 + 1.0 / s; (void) x1;
    LG = R::lgammafn(x1);
    t1 = std::exp(LG*s);
  };
  const bool scalar = n_mean == 1 && n_sigma == 1;
  if (scalar) set_params(0);
  for (int i = 0; i < n; i++) {
    if (!scalar) set_params(i);
    const double yy = y[i];
    const double C = std::log(yy / m);
    const double t0 = std::exp(C*s);
    const double t2 = t0*t1;
    o_y[i] = (-std::pow(s, 4)*t2 + 6*std::pow(s, 3)*t0*t1 - 11*std::pow(s, 2)*t2 + 6*s*t0*t1 - 6*s + 6)/std::pow(yy, 4);
  }
  (void) LG;
  return List::create(Named("y") = o_y);
}
