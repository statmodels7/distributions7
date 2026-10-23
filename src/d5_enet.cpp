#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "pt_enet_logz.h"
using namespace Rcpp;

// The fifth derivatives of the elastic-net log-density in (mu, lambda, alpha):
//   l = -a s d - c d^2/2 - log Z(a, c),  a = lambda alpha, c = lambda (1 - alpha),
//   d = y - mu, s = sign(d), Z the normalizing constant of pt_enet_logz.h.
// log Z is a function of (a, c) alone, and d^(i+j) log Z / da^i dc^j for
// i + j <= 5 comes from d5_enet_logz5() below, which is enet_logz_ac() of
// pt_enet_logz.h run to order five (that header stops at four and its
// array is sized for it). The map (lambda, alpha) -> (a, c) is bilinear, so
// the chain onto the reported parameters is written out per component. The
// loop body is derived in misure_gas_congiunto/c/fase4/D5/enet/derive.py.

namespace {

// d[6 i + j] = d^(i+j) log Z / da^i dc^j for 1 <= i + j <= 5, as
// enet_logz_ac() in pt_enet_logz.h, with the same two representations
// (x < 2 through L = log M0, x >= 2 through t = c/a^2 and F(t)), to order five.
inline void d5_enet_logz5(double a, double c, double* d) {
  const int n = 5;
  const double x = a / std::sqrt(c);
  for (int k = 0; k < 36; ++k) d[k] = 0.0;
  if (x >= 2.0) {
    const double t = 1.0 / (x * x);
    double m[12], mu[6], F[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    d7::enet_s_moments(x, 2 * n, m);
    mu[0] = 1.0;
    double h = 1.0, x2 = 1.0;
    for (int j = 1; j <= n; ++j) {
      h *= -0.5;
      x2 *= x * x;
      mu[j] = h * x2 * m[2 * j];
    }
    d7::enet_cumulants(mu, n, F);
    for (int ord = 1; ord <= n; ++ord) {
      for (int j = 0; j <= ord; ++j) {
        const int i = ord - j;
        double b[8] = {1, 0, 0, 0, 0, 0, 0, 0};
        int len = 1;
        for (int q = 0; q < i; ++q) len = d7::enet_euler_step(b, len, j, 2.0, q);
        double v = 0.0, tm = 1.0;
        for (int mm = 0; mm < len; ++mm) {
          if (b[mm] != 0.0) v += b[mm] * tm * F[j + mm];
          tm *= t;
        }
        v *= std::pow(a, -2.0 * j - i);
        if (j == 0) {
          double f = 1.0;
          for (int k = 2; k < i; ++k) f *= k;
          v += ((i % 2) ? -f : f) * std::pow(a, -i);
        }
        d[6 * i + j] = v;
      }
    }
    return;
  }
  double L[6];
  d7::enet_logM_derivs(x, n, L);
  for (int ord = 1; ord <= n; ++ord) {
    for (int j = 0; j <= ord; ++j) {
      const int i = ord - j;
      double b[8] = {1, 0, 0, 0, 0, 0, 0, 0};
      int len = 1;
      for (int q = 0; q < j; ++q) len = d7::enet_euler_step(b, len, i, 0.5, q);
      double v = 0.0, xm = 1.0;
      for (int mm = 0; mm < len; ++mm) {
        if (i + mm >= 1 && b[mm] != 0.0) v += b[mm] * xm * L[i + mm];
        xm *= x;
      }
      v *= std::pow(c, -0.5 * i - j);
      if (i == 0) {
        double f = 1.0;
        for (int k = 2; k < j; ++k) f *= k;
        v += -0.5 * ((j % 2) ? f : -f) * std::pow(c, -j);
      }
      d[6 * i + j] = v;
    }
  }
}

}  // namespace

// [[Rcpp::export]]
List enet_deriv5_cpp(NumericVector y, NumericVector mu, NumericVector lambda,
                     NumericVector alpha) {
  const int n = std::max(y.size(), std::max(mu.size(), std::max(lambda.size(), alpha.size())));
  const int n_y = y.size(), n_mu = mu.size(), n_l = lambda.size(), n_a = alpha.size();
  NumericVector o_mu_mu_mu_mu_mu(n), o_mu_mu_mu_mu_lambda(n), o_mu_mu_mu_mu_alpha(n),
    o_mu_mu_mu_lambda_lambda(n), o_mu_mu_mu_lambda_alpha(n), o_mu_mu_mu_alpha_alpha(n),
    o_mu_mu_lambda_lambda_lambda(n), o_mu_mu_lambda_lambda_alpha(n),
    o_mu_mu_lambda_alpha_alpha(n), o_mu_mu_alpha_alpha_alpha(n),
    o_mu_lambda_lambda_lambda_lambda(n), o_mu_lambda_lambda_lambda_alpha(n),
    o_mu_lambda_lambda_alpha_alpha(n), o_mu_lambda_alpha_alpha_alpha(n),
    o_mu_alpha_alpha_alpha_alpha(n), o_lambda_lambda_lambda_lambda_lambda(n),
    o_lambda_lambda_lambda_lambda_alpha(n), o_lambda_lambda_lambda_alpha_alpha(n),
    o_lambda_lambda_alpha_alpha_alpha(n), o_lambda_alpha_alpha_alpha_alpha(n),
    o_alpha_alpha_alpha_alpha_alpha(n);
  for (int i = 0; i < n; i++) {
    const double yy = y[i % n_y];
    const double mm = mu[i % n_mu];
    const double lam = lambda[i % n_l];
    const double al = alpha[i % n_a];
    const double aa = lam*al;
    const double cc = lam*(1.0 - al);
    double Zd[36];
    d5_enet_logz5(aa, cc, Zd);
    const double Z_0_3 = Zd[3];
    const double Z_0_4 = Zd[4];
    const double Z_0_5 = Zd[5];
    const double Z_1_2 = Zd[8];
    const double Z_1_3 = Zd[9];
    const double Z_1_4 = Zd[10];
    const double Z_2_1 = Zd[13];
    const double Z_2_2 = Zd[14];
    const double Z_2_3 = Zd[15];
    const double Z_3_0 = Zd[18];
    const double Z_3_1 = Zd[19];
    const double Z_3_2 = Zd[20];
    const double Z_4_0 = Zd[24];
    const double Z_4_1 = Zd[25];
    const double Z_5_0 = Zd[30];
        const double t0 = std::pow(al, 4);
        const double t1 = al - 1;
        const double t2 = -t1;
        const double t3 = t1*std::pow(t2, 3);
        const double t4 = std::pow(al, 3);
        const double t5 = std::pow(al, 2);
        const double t6 = al*t2;
        const double t7 = al*t1 - t6;
        const double t8 = al*t7 - t2*t5;
        const double t9 = al*t8 - t2*t4;
        const double t10 = std::pow(t2, 2);
        const double t11 = al*t1;
        const double t12 = t11*t2 + t2*t7;
        const double t13 = t10*t11 + t12*t2;
        const double t14 = al*t12 + t2*t8;
        const double t15 = Z_4_1*lam;
        const double t16 = Z_5_0*lam - t15;
        const double t17 = Z_0_4*t2 + Z_1_3*al;
        const double t18 = Z_1_4*lam;
        const double t19 = -Z_0_5*lam + t18;
        const double t20 = 2*al;
        const double t21 = t20 - 2;
        const double t22 = t1*t21;
        const double t23 = 4*al;
        const double t24 = t23 - 2;
        const double t25 = al*t24 + t11 + t5 - 3*t6;
        const double t26 = t1*t2;
        const double t27 = -2*t11 + t2*t24 + t26 + 2*t6;
        const double t28 = Z_3_2*lam;
        const double t29 = t15 - t28;
        const double t30 = Z_2_3*lam;
        const double t31 = -t18 + t30;
        const double t32 = t28 - t30;
        const double t33 = 2*t2;
        const double t34 = -Z_0_3 + Z_1_2;
        const double t35 = Z_3_1*lam;
        const double t36 = Z_4_0*lam - t35;
        const double t37 = 6*t36;
        const double t38 = lam*t4;
        const double t39 = Z_1_3*lam;
        const double t40 = -Z_0_4*lam + t39;
        const double t41 = lam*t1;
        const double t42 = t10*t41;
        const double t43 = Z_2_2*lam;
        const double t44 = t35 - t43;
        const double t45 = lam*t8;
        const double t46 = -t39 + t43;
        const double t47 = -lam*t12;
        const double t48 = 6*lam;
        const double t49 = 18*lam;
        const double t50 = al*lam;
        const double t51 = std::pow(lam, 2);
        const double t52 = t5*t51;
        const double t53 = lam*t2;
        const double t54 = t26*t51;
        const double t55 = lam*t24;
        const double t56 = 2*t55;
        const double t57 = -lam*t21 - t41 + t53;
        const double t58 = lam*t7;
        const double t59 = lam*(lam*t5 + t58);
        const double t60 = lam*(lam*t26 - t58);
        const double t61 = 4*t51;
        const double t62 = 12*t51;
        const double t63 = std::pow(lam, 3);
        const double t64 = al*t63;
        const double t65 = t1*t63;
        const double t66 = lam*(t41 + t50);
        const double t67 = lam*(al*t51 + t66);
        const double t68 = -lam*(-t1*t51 - t66);
        const double t69 = std::pow(lam, 4);
        o_mu_mu_mu_mu_mu[i] = 0;
        o_mu_mu_mu_mu_lambda[i] = 0;
        o_mu_mu_mu_mu_alpha[i] = 0;
        o_mu_mu_mu_lambda_lambda[i] = 0;
        o_mu_mu_mu_lambda_alpha[i] = 0;
        o_mu_mu_mu_alpha_alpha[i] = 0;
        o_mu_mu_lambda_lambda_lambda[i] = 0;
        o_mu_mu_lambda_lambda_alpha[i] = 0;
        o_mu_mu_lambda_alpha_alpha[i] = 0;
        o_mu_mu_alpha_alpha_alpha[i] = 0;
        o_mu_lambda_lambda_lambda_lambda[i] = 0;
        o_mu_lambda_lambda_lambda_alpha[i] = 0;
        o_mu_lambda_lambda_alpha_alpha[i] = 0;
        o_mu_lambda_alpha_alpha_alpha[i] = 0;
        o_mu_alpha_alpha_alpha_alpha[i] = 0;
        o_lambda_lambda_lambda_lambda_lambda[i] = -t0*(Z_4_1*t2 + Z_5_0*al) + t13*(Z_1_4*t2 + Z_2_3*al) + t14*(Z_2_3*t2 + Z_3_2*al) + t3*(Z_0_5*t2 + Z_1_4*al) + t9*(Z_3_2*t2 + Z_4_1*al);
        o_lambda_lambda_lambda_lambda_alpha[i] = -t0*t16 + t1*t10*(-Z_0_4 + Z_1_3) + t10*t17 + t12*(-Z_1_3 + Z_2_2) + t13*t31 + t14*t32 + t17*t22 + t19*t3 + t25*(Z_2_2*t2 + Z_3_1*al) + t27*(Z_1_3*t2 + Z_2_2*al) + t29*t9 - t4*(-Z_3_1 + Z_4_0) - 3*t5*(Z_3_1*t2 + Z_4_0*al) + t8*(-Z_2_2 + Z_3_1);
        o_lambda_lambda_lambda_alpha_alpha[i] = -Z_0_3*t33 + 4*Z_1_2*t2 - Z_1_2*t20 + Z_2_1*t23 - Z_2_1*t33 - Z_3_0*t20 - 2*t1*t34 + t10*t40 - t16*t38 - t19*t42 + t22*t40 - t23*(-Z_2_1 + Z_3_0) + 2*t24*(-Z_1_2 + Z_2_1) + 2*t25*t44 + 2*t27*t46 + t29*(t38 + t45) + t31*(t42 + t47) + t32*(-t45 - t47) + t33*t34 - t37*t5 + t40*(t10 - 2*t26);
        o_lambda_lambda_alpha_alpha_alpha[i] = Z_0_3*t48 - Z_1_2*t49 + Z_2_1*t49 - Z_3_0*t48 - t16*t52 + t19*t54 + t29*(t52 + t59) + t31*(-t54 - t60) + t32*(-t59 + t60) - t37*t50 + t40*t41 - t40*t53 - t40*t57 + t44*(lam*t20 + t55) + t44*(lam*t23 + t56) + t46*(-t56 + t57) + t46*(-t41 + t53 - t55);
        o_lambda_alpha_alpha_alpha_alpha[i] = -t16*t64 - t19*t65 + t29*(t64 + t67) + t31*(t65 + t68) + t32*(-t67 - t68) - t36*t61 + t40*t61 + t44*t62 - t46*t62;
        o_alpha_alpha_alpha_alpha_alpha[i] = -t16*t69 - t19*t69 + 4*t29*t69 + 4*t31*t69 - 6*t32*t69;
  }
  return List::create(Named("mu_mu_mu_mu_mu") = o_mu_mu_mu_mu_mu,
                      Named("mu_mu_mu_mu_lambda") = o_mu_mu_mu_mu_lambda,
                      Named("mu_mu_mu_mu_alpha") = o_mu_mu_mu_mu_alpha,
                      Named("mu_mu_mu_lambda_lambda") = o_mu_mu_mu_lambda_lambda,
                      Named("mu_mu_mu_lambda_alpha") = o_mu_mu_mu_lambda_alpha,
                      Named("mu_mu_mu_alpha_alpha") = o_mu_mu_mu_alpha_alpha,
                      Named("mu_mu_lambda_lambda_lambda") = o_mu_mu_lambda_lambda_lambda,
                      Named("mu_mu_lambda_lambda_alpha") = o_mu_mu_lambda_lambda_alpha,
                      Named("mu_mu_lambda_alpha_alpha") = o_mu_mu_lambda_alpha_alpha,
                      Named("mu_mu_alpha_alpha_alpha") = o_mu_mu_alpha_alpha_alpha,
                      Named("mu_lambda_lambda_lambda_lambda") = o_mu_lambda_lambda_lambda_lambda,
                      Named("mu_lambda_lambda_lambda_alpha") = o_mu_lambda_lambda_lambda_alpha,
                      Named("mu_lambda_lambda_alpha_alpha") = o_mu_lambda_lambda_alpha_alpha,
                      Named("mu_lambda_alpha_alpha_alpha") = o_mu_lambda_alpha_alpha_alpha,
                      Named("mu_alpha_alpha_alpha_alpha") = o_mu_alpha_alpha_alpha_alpha,
                      Named("lambda_lambda_lambda_lambda_lambda") = o_lambda_lambda_lambda_lambda_lambda,
                      Named("lambda_lambda_lambda_lambda_alpha") = o_lambda_lambda_lambda_lambda_alpha,
                      Named("lambda_lambda_lambda_alpha_alpha") = o_lambda_lambda_lambda_alpha_alpha,
                      Named("lambda_lambda_alpha_alpha_alpha") = o_lambda_lambda_alpha_alpha_alpha,
                      Named("lambda_alpha_alpha_alpha_alpha") = o_lambda_alpha_alpha_alpha_alpha,
                      Named("alpha_alpha_alpha_alpha_alpha") = o_alpha_alpha_alpha_alpha_alpha);
}
