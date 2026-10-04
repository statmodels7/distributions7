#ifndef D7_PT_LOGNORMAL2_H
#define D7_PT_LOGNORMAL2_H

#include <Rcpp.h>
#include <cmath>

// The lognormal in the mean m and the variance v of Y: one function per
// component and order for the quantities the scalar registry reads, lifted
// from the generated kernels in lognormal2.cpp (which records the
// derivation) and called by them and by the registry (d7_ccallable.cpp).
// The caller passes S = log1p(v/m^2) and C = log(y/m).

namespace d7 {

inline double lognormal2_score_mean(double m, double v, double S, double C) {
  const double t0 = 1.0/m;
  const double t1 = 1.0/S;
  const double t2 = t0*t1;
  const double t3 = 1.0/(std::pow(m, 2) + v);
  const double t4 = m*t3;
  const double t6 = std::pow(S, -2);
  const double t5 = std::pow(C, 2);
  const double t7 = t5*t6;
  return C*t2 - t0*t7 + (3.0/4.0)*t0 - t1*t4 + t2 + t4*t7 - 1.0/4.0*t4;
}

inline double lognormal2_score_var(double m, double v, double S, double C) {
  const double t1 = 1.0/S;
  const double t3 = 1.0/(std::pow(m, 2) + v);
  const double t6 = std::pow(S, -2);
  const double t5 = std::pow(C, 2);
  return (1.0/8.0)*t3*(-4*t1 + 4*t5*t6 - 1);
}

inline double lognormal2_hess_mean_mean(double m, double v, double S, double C) {
  const double t0 = std::pow(m, 2);
  const double t1 = 1.0/t0;
  const double t2 = std::pow(S, -2);
  const double t3 = 1.0/S;
  const double t4 = 1.0/(t0 + v);
  const double t6 = std::pow(S, -3);
  const double t9 = t2*t4;
  const double t10 = 4*t9;
  const double t11 = 2*t0;
  const double t12 = 1.0/(std::pow(m, 4) + t11*v + std::pow(v, 2));
  const double t5 = std::pow(C, 2);
  const double t7 = t5*t6;
  const double t8 = 4*t7;
  return 4*C*t1*t2 - C*t1*t3 - C*t10 + 2*t0*t12*t2 + 2*t0*t12*t3 - t0*t12*t8 + (1.0/2.0)*t0*t12 + t1*t2*t5 + 2*t1*t2 - 2*t1*t3 - t1*t8 - 3.0/4.0*t1 - t10 - t11*t12*t2*t5 + t2*t4*t5 - t3*t4 + 8*t4*t5*t6 - 1.0/4.0*t4;
}

inline double lognormal2_hess_var_var(double m, double v, double S, double C) {
  const double t0 = std::pow(m, 2);
  const double t2 = std::pow(S, -2);
  const double t3 = 1.0/S;
  const double t6 = std::pow(S, -3);
  const double t11 = 2*t0;
  const double t12 = 1.0/(std::pow(m, 4) + t11*v + std::pow(v, 2));
  const double t13 = (1.0/2.0)*t2;
  const double t5 = std::pow(C, 2);
  const double t7 = t5*t6;
  return t12*(-t13*t5 + t13 + (1.0/2.0)*t3 - t7 + 1.0/8.0);
}

inline double lognormal2_expected_mean_mean(double m, double v, double S) {
  const double t0 = std::pow(m, 2);
  const double t1 = 1.0/(t0 + v);
  const double t2 = 1.0/t0;
  const double t3 = 1.0/S;
  const double t4 = 2*t3;
  const double t5 = std::pow(v, 2);
  const double t6 = 1.0/(std::pow(m, 4) + 2*t0*v + t5);
  const double t7 = t6*v;
  const double t8 = t0*t6;
  const double t9 = t2*t5*t6;
  const double t10 = std::pow(S, -2);
  return -t1*t2*t4*v + (1.0/2.0)*t1*t2*v + (1.0/2.0)*t1 - 2*t10*t9 - t3*t8 - t4*t7 - t4*t9 - t7 - 1.0/2.0*t8 - 1.0/2.0*t9;
}

inline double lognormal2_expected_var_var(double m, double v, double S) {
  const double t0 = std::pow(m, 2);
  const double t3 = 1.0/S;
  const double t4 = 2*t3;
  const double t5 = std::pow(v, 2);
  const double t6 = 1.0/(std::pow(m, 4) + 2*t0*v + t5);
  const double t10 = std::pow(S, -2);
  return (1.0/8.0)*t6*(-4*t10 - t4);
}

inline double lognormal2_dexpected_mean_mean_mean(double m, double v, double S) {
  const double t0 = 1.0/m;
  const double t1 = std::pow(m, 2);
  const double t2 = 1.0/(t1 + v);
  const double t3 = t0*t2;
  const double t4 = 1.0/S;
  const double t5 = 2*t4;
  const double t6 = std::pow(S, -2);
  const double t7 = 4*v;
  const double t8 = std::pow(m, 3);
  const double t9 = 1.0/t8;
  const double t10 = t2*t9;
  const double t11 = std::pow(m, 4);
  const double t12 = std::pow(v, 2);
  const double t13 = 2*t1;
  const double t14 = 1.0/(t11 + t12 + t13*v);
  const double t15 = m*t14;
  const double t16 = 2*t6;
  const double t17 = t14*t4;
  const double t18 = 2*m;
  const double t19 = 5*v;
  const double t20 = t0*t14;
  const double t21 = t12*t9;
  const double t22 = t14*t21;
  const double t23 = t14*t6;
  const double t24 = t0*t23;
  const double t25 = t0*t17;
  const double t26 = std::pow(S, -3);
  const double t27 = 8*t26;
  const double t28 = 4*t6;
  const double t29 = std::pow(m, 6);
  const double t30 = std::pow(v, 3);
  const double t31 = t11*v;
  const double t32 = t1*t12;
  const double t33 = 1.0/(t29 + t30 + 3*t31 + 3*t32);
  const double t34 = t33*t6;
  const double t35 = m*v;
  const double t36 = t33*t4;
  const double t37 = t36*v;
  const double t38 = 2*t34;
  const double t39 = 4*t8;
  const double t40 = t26*t33;
  const double t41 = t0*t12;
  const double t42 = std::pow(m, 5);
  const double t43 = std::pow(m, 8);
  const double t44 = std::pow(v, 4);
  const double t45 = 4*t1;
  const double t46 = 1.0/(6*t11*t12 + t29*t7 + t30*t45 + t43 + t44);
  const double t47 = 3*t46;
  const double t48 = t12*t46;
  const double t49 = m*t48;
  const double t50 = 10*t8;
  const double t51 = t46*v;
  const double t52 = t44*t46*t9;
  const double t53 = t46*t6;
  const double t54 = 4*t42;
  const double t55 = t4*t46;
  const double t56 = t30*t46;
  const double t57 = t0*t56;
  const double t58 = 16*t26;
  const double t59 = t12*t53;
  const double t60 = 4*t59;
  const double t61 = t53*v;
  const double t62 = t12*t55;
  const double t63 = t55*t8;
  const double t64 = std::pow(m, 7);
  const double t65 = 1.0/(std::pow(m, 10) + 5*t1*t44 + 10*t11*t30 + 10*t12*t29 + t19*t43 + std::pow(v, 5));
  const double t66 = t6*t65;
  const double t67 = 4*t66;
  const double t68 = t4*t65;
  const double t69 = t26*t65;
  const double t70 = t30*t69;
  const double t71 = t30*t66;
  const double t72 = t42*t66;
  const double t73 = t30*t68;
  const double t74 = t42*t68;
  const double t75 = t0*t44;
  const double t76 = 16*t69;
  const double t77 = t12*t8;
  const double t78 = 8*t66;
  return -13.0/2.0*m*t37 - m*t60 + (43.0/2.0)*m*t62 + 32*m*t70 + 18*m*t71 + (11.0/2.0)*m*t73 - t10*t5*v + t10*t6*t7 + t15*t16 - 4*t15 - t16*t57 + t17*t18 + (15.0/2.0)*t17*t21 - t19*t20 + t22*t27 - t22*t28 - t22 - 14*t24*v + (25.0/2.0)*t25*v - 32*t26*t57 - t28*t52 - t3*t5 + t3 + 10*t34*t35 + 4*t34*t41 - t36*t39 - 7.0/2.0*t36*t41 - t38*t8 + (5.0/2.0)*t4*t52 + (29.0/2.0)*t4*t57 - 8*t40*t41 + t42*t47 - t49*t58 + 12*t49 + t50*t51 - t50*t61 - t52*t58 + t52 - t53*t54 + t54*t55 + 6*t57 + (27.0/2.0)*t63*v + t64*t67 + 2*t64*t68 + 16*t66*t77 + (3.0/2.0)*t68*t75 + (17.0/2.0)*t68*t77 + 10*t72*v + (13.0/2.0)*t74*v + t75*t76 + t75*t78 + t76*t77;
}

inline double lognormal2_dexpected_var_var_var(double m, double v, double S) {
  const double t1 = std::pow(m, 2);
  const double t4 = 1.0/S;
  const double t6 = std::pow(S, -2);
  const double t7 = 4*v;
  const double t11 = std::pow(m, 4);
  const double t12 = std::pow(v, 2);
  const double t19 = 5*v;
  const double t26 = std::pow(S, -3);
  const double t29 = std::pow(m, 6);
  const double t30 = std::pow(v, 3);
  const double t31 = t11*v;
  const double t32 = t1*t12;
  const double t33 = 1.0/(t29 + t30 + 3*t31 + 3*t32);
  const double t34 = t33*t6;
  const double t36 = t33*t4;
  const double t40 = t26*t33;
  const double t43 = std::pow(m, 8);
  const double t44 = std::pow(v, 4);
  const double t45 = 4*t1;
  const double t46 = 1.0/(6*t11*t12 + t29*t7 + t30*t45 + t43 + t44);
  const double t53 = t46*t6;
  const double t55 = t4*t46;
  const double t61 = t53*v;
  const double t65 = 1.0/(std::pow(m, 10) + 5*t1*t44 + 10*t11*t30 + 10*t12*t29 + t19*t43 + std::pow(v, 5));
  const double t66 = t6*t65;
  const double t68 = t4*t65;
  const double t69 = t26*t65;
  const double t86 = t1*v;
  const double t94 = (1.0/2.0)*t55;
  const double t95 = t94*v;
  const double t100 = (5.0/2.0)*t66;
  const double t101 = 2*t69;
  const double t102 = (5.0/4.0)*t66;
  const double t103 = (1.0/8.0)*t68;
  const double t104 = t69*t7;
  return t1*t104 + t1*t53 + t1*t94 + t100*t86 + t101*t11 + t101*t12 + t102*t11 + t102*t12 + t103*t11 + t103*t12 - t34 - 1.0/8.0*t36 - t40 + t61 + (1.0/4.0)*t68*t86 + t95;
}

inline void lognormal2_score_curv(int k, double y, const double* th,
                                  double* out) {
  const double m = th[0], v = th[1];
  const double S = std::log1p(v / (m * m));
  const double C = std::log(y / m);
  if (k == 0) {
    out[0] = lognormal2_score_mean(m, v, S, C);
    out[1] = lognormal2_hess_mean_mean(m, v, S, C);
  } else {
    out[0] = lognormal2_score_var(m, v, S, C);
    out[1] = lognormal2_hess_var_var(m, v, S, C);
  }
}

inline void lognormal2_info_dinfo(int k, double y, const double* th,
                                  double* out) {
  const double m = th[0], v = th[1];
  const double S = std::log1p(v / (m * m));
  if (k == 0) {
    out[0] = lognormal2_expected_mean_mean(m, v, S);
    out[1] = lognormal2_dexpected_mean_mean_mean(m, v, S);
  } else {
    out[0] = lognormal2_expected_var_var(m, v, S);
    out[1] = lognormal2_dexpected_var_var_var(m, v, S);
  }
}

} // namespace d7

#endif
