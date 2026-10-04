#ifndef D7_PT_LOGNORMAL2_H
#define D7_PT_LOGNORMAL2_H

#include <Rcpp.h>
#include <cmath>

// The lognormal in the mean m and the variance v of Y: one function per
// component and order for the quantities the scalar registry reads, lifted
// from the generated kernels in lognormal2.cpp (which records the
// derivation) and called by them and by the registry (d7_ccallable.cpp).
// The quantities of the parameters alone, S = log1p(v/m^2) among them, are
// formed once per parameter value by lognormal2_<kernel>_par(); the caller
// passes that struct and C = log(y/m).

namespace d7 {

struct Lognormal2GradientPar { double S, t0, t1, t2, t3, t4, t6; };

inline Lognormal2GradientPar lognormal2_gradient_par(double m, double v) {
  const double S = std::log1p(v / (m * m));
  const double t0 = 1.0/m;
  const double t1 = 1.0/S;
  const double t2 = t0*t1;
  const double t3 = 1.0/(std::pow(m, 2) + v);
  const double t4 = m*t3;
  const double t6 = std::pow(S, -2);
  Lognormal2GradientPar P;
  P.S = S;
  P.t0 = t0;
  P.t1 = t1;
  P.t2 = t2;
  P.t3 = t3;
  P.t4 = t4;
  P.t6 = t6;
  return P;
}

struct Lognormal2HessianPar { double S, t0, t1, t2, t3, t4, t6, t9, t10, t11, t12, t13, t14, t15, t16; };

inline Lognormal2HessianPar lognormal2_hessian_par(double m, double v) {
  const double S = std::log1p(v / (m * m));
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
  const double t13 = (1.0/2.0)*t2;
  const double t14 = 1.0/m;
  const double t15 = t14*t9;
  const double t16 = m*t12;
  Lognormal2HessianPar P;
  P.S = S;
  P.t0 = t0;
  P.t1 = t1;
  P.t2 = t2;
  P.t3 = t3;
  P.t4 = t4;
  P.t6 = t6;
  P.t9 = t9;
  P.t10 = t10;
  P.t11 = t11;
  P.t12 = t12;
  P.t13 = t13;
  P.t14 = t14;
  P.t15 = t15;
  P.t16 = t16;
  return P;
}

struct Lognormal2ExpectedPar { double S, t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12; };

inline Lognormal2ExpectedPar lognormal2_expected_par(double m, double v) {
  const double S = std::log1p(v / (m * m));
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
  const double t11 = 1.0/m;
  const double t12 = t11*t3;
  Lognormal2ExpectedPar P;
  P.S = S;
  P.t0 = t0;
  P.t1 = t1;
  P.t2 = t2;
  P.t3 = t3;
  P.t4 = t4;
  P.t5 = t5;
  P.t6 = t6;
  P.t7 = t7;
  P.t8 = t8;
  P.t9 = t9;
  P.t10 = t10;
  P.t11 = t11;
  P.t12 = t12;
  return P;
}

struct Lognormal2DexpectedPar { double S, t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t53, t54, t55, t56, t57, t58, t59, t60, t61, t62, t63, t64, t65, t66, t67, t68, t69, t70, t71, t72, t73, t74, t75, t76, t77, t78, t79, t80, t81, t82, t83, t84, t85, t86, t87, t88, t89, t90, t91, t92, t93, t94, t95, t96, t97, t98, t99, t100, t101, t102, t103, t104, t105, t106, t107; };

inline Lognormal2DexpectedPar lognormal2_dexpected_par(double m, double v) {
  const double S = std::log1p(v / (m * m));
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
  const double t79 = 1.0/t1;
  const double t80 = t2*t79;
  const double t81 = (3.0/2.0)*t14;
  const double t82 = t79*v;
  const double t83 = t17*t82;
  const double t84 = t12*t79;
  const double t85 = t11*t46;
  const double t86 = t1*v;
  const double t87 = t56*t79;
  const double t88 = t55*t86;
  const double t89 = t44*t79;
  const double t90 = 8*t69;
  const double t91 = m*t38;
  const double t92 = (1.0/4.0)*m;
  const double t93 = t36*t92;
  const double t94 = (1.0/2.0)*t55;
  const double t95 = t94*v;
  const double t96 = (1.0/2.0)*t8;
  const double t97 = 4*t0;
  const double t98 = t0*t62;
  const double t99 = m*t12;
  const double t100 = (5.0/2.0)*t66;
  const double t101 = 2*t69;
  const double t102 = (5.0/4.0)*t66;
  const double t103 = (1.0/8.0)*t68;
  const double t104 = t69*t7;
  const double t105 = 2*t59;
  const double t106 = (1.0/2.0)*t4;
  const double t107 = t12*t66;
  Lognormal2DexpectedPar P;
  P.S = S;
  P.t0 = t0;
  P.t1 = t1;
  P.t2 = t2;
  P.t3 = t3;
  P.t4 = t4;
  P.t5 = t5;
  P.t6 = t6;
  P.t7 = t7;
  P.t8 = t8;
  P.t9 = t9;
  P.t10 = t10;
  P.t11 = t11;
  P.t12 = t12;
  P.t13 = t13;
  P.t14 = t14;
  P.t15 = t15;
  P.t16 = t16;
  P.t17 = t17;
  P.t18 = t18;
  P.t19 = t19;
  P.t20 = t20;
  P.t21 = t21;
  P.t22 = t22;
  P.t23 = t23;
  P.t24 = t24;
  P.t25 = t25;
  P.t26 = t26;
  P.t27 = t27;
  P.t28 = t28;
  P.t29 = t29;
  P.t30 = t30;
  P.t31 = t31;
  P.t32 = t32;
  P.t33 = t33;
  P.t34 = t34;
  P.t35 = t35;
  P.t36 = t36;
  P.t37 = t37;
  P.t38 = t38;
  P.t39 = t39;
  P.t40 = t40;
  P.t41 = t41;
  P.t42 = t42;
  P.t43 = t43;
  P.t44 = t44;
  P.t45 = t45;
  P.t46 = t46;
  P.t47 = t47;
  P.t48 = t48;
  P.t49 = t49;
  P.t50 = t50;
  P.t51 = t51;
  P.t52 = t52;
  P.t53 = t53;
  P.t54 = t54;
  P.t55 = t55;
  P.t56 = t56;
  P.t57 = t57;
  P.t58 = t58;
  P.t59 = t59;
  P.t60 = t60;
  P.t61 = t61;
  P.t62 = t62;
  P.t63 = t63;
  P.t64 = t64;
  P.t65 = t65;
  P.t66 = t66;
  P.t67 = t67;
  P.t68 = t68;
  P.t69 = t69;
  P.t70 = t70;
  P.t71 = t71;
  P.t72 = t72;
  P.t73 = t73;
  P.t74 = t74;
  P.t75 = t75;
  P.t76 = t76;
  P.t77 = t77;
  P.t78 = t78;
  P.t79 = t79;
  P.t80 = t80;
  P.t81 = t81;
  P.t82 = t82;
  P.t83 = t83;
  P.t84 = t84;
  P.t85 = t85;
  P.t86 = t86;
  P.t87 = t87;
  P.t88 = t88;
  P.t89 = t89;
  P.t90 = t90;
  P.t91 = t91;
  P.t92 = t92;
  P.t93 = t93;
  P.t94 = t94;
  P.t95 = t95;
  P.t96 = t96;
  P.t97 = t97;
  P.t98 = t98;
  P.t99 = t99;
  P.t100 = t100;
  P.t101 = t101;
  P.t102 = t102;
  P.t103 = t103;
  P.t104 = t104;
  P.t105 = t105;
  P.t106 = t106;
  P.t107 = t107;
  return P;
}

inline double lognormal2_score_mean(const Lognormal2GradientPar& P, double C) {
  const double t0 = P.t0;
  const double t1 = P.t1;
  const double t2 = P.t2;
  const double t4 = P.t4;
  const double t6 = P.t6;
  const double t5 = std::pow(C, 2);
  const double t7 = t5*t6;
  return C*t2 - t0*t7 + (3.0/4.0)*t0 - t1*t4 + t2 + t4*t7 - 1.0/4.0*t4;
}

inline double lognormal2_score_var(const Lognormal2GradientPar& P, double C) {
  const double t1 = P.t1;
  const double t3 = P.t3;
  const double t6 = P.t6;
  const double t5 = std::pow(C, 2);
  return (1.0/8.0)*t3*(-4*t1 + 4*t5*t6 - 1);
}

inline double lognormal2_hess_mean_mean(const Lognormal2HessianPar& P, double C) {
  const double t0 = P.t0;
  const double t1 = P.t1;
  const double t2 = P.t2;
  const double t3 = P.t3;
  const double t4 = P.t4;
  const double t6 = P.t6;
  const double t10 = P.t10;
  const double t11 = P.t11;
  const double t12 = P.t12;
  const double t5 = std::pow(C, 2);
  const double t7 = t5*t6;
  const double t8 = 4*t7;
  return 4*C*t1*t2 - C*t1*t3 - C*t10 + 2*t0*t12*t2 + 2*t0*t12*t3 - t0*t12*t8 + (1.0/2.0)*t0*t12 + t1*t2*t5 + 2*t1*t2 - 2*t1*t3 - t1*t8 - 3.0/4.0*t1 - t10 - t11*t12*t2*t5 + t2*t4*t5 - t3*t4 + 8*t4*t5*t6 - 1.0/4.0*t4;
}

inline double lognormal2_hess_var_var(const Lognormal2HessianPar& P, double C) {
  const double t3 = P.t3;
  const double t6 = P.t6;
  const double t12 = P.t12;
  const double t13 = P.t13;
  const double t5 = std::pow(C, 2);
  const double t7 = t5*t6;
  return t12*(-t13*t5 + t13 + (1.0/2.0)*t3 - t7 + 1.0/8.0);
}

inline double lognormal2_expected_mean_mean(const Lognormal2ExpectedPar& P, double v) {
  const double t1 = P.t1;
  const double t2 = P.t2;
  const double t3 = P.t3;
  const double t4 = P.t4;
  const double t7 = P.t7;
  const double t8 = P.t8;
  const double t9 = P.t9;
  const double t10 = P.t10;
  return -t1*t2*t4*v + (1.0/2.0)*t1*t2*v + (1.0/2.0)*t1 - 2*t10*t9 - t3*t8 - t4*t7 - t4*t9 - t7 - 1.0/2.0*t8 - 1.0/2.0*t9;
}

inline double lognormal2_expected_var_var(const Lognormal2ExpectedPar& P) {
  const double t4 = P.t4;
  const double t6 = P.t6;
  const double t10 = P.t10;
  return (1.0/8.0)*t6*(-4*t10 - t4);
}

inline double lognormal2_dexpected_mean_mean_mean(const Lognormal2DexpectedPar& P, double m, double v) {
  const double t3 = P.t3;
  const double t4 = P.t4;
  const double t5 = P.t5;
  const double t6 = P.t6;
  const double t7 = P.t7;
  const double t8 = P.t8;
  const double t10 = P.t10;
  const double t15 = P.t15;
  const double t16 = P.t16;
  const double t17 = P.t17;
  const double t18 = P.t18;
  const double t19 = P.t19;
  const double t20 = P.t20;
  const double t21 = P.t21;
  const double t22 = P.t22;
  const double t24 = P.t24;
  const double t25 = P.t25;
  const double t26 = P.t26;
  const double t27 = P.t27;
  const double t28 = P.t28;
  const double t34 = P.t34;
  const double t35 = P.t35;
  const double t36 = P.t36;
  const double t37 = P.t37;
  const double t38 = P.t38;
  const double t39 = P.t39;
  const double t40 = P.t40;
  const double t41 = P.t41;
  const double t42 = P.t42;
  const double t47 = P.t47;
  const double t49 = P.t49;
  const double t50 = P.t50;
  const double t51 = P.t51;
  const double t52 = P.t52;
  const double t53 = P.t53;
  const double t54 = P.t54;
  const double t55 = P.t55;
  const double t57 = P.t57;
  const double t58 = P.t58;
  const double t60 = P.t60;
  const double t61 = P.t61;
  const double t62 = P.t62;
  const double t63 = P.t63;
  const double t64 = P.t64;
  const double t66 = P.t66;
  const double t67 = P.t67;
  const double t68 = P.t68;
  const double t70 = P.t70;
  const double t71 = P.t71;
  const double t72 = P.t72;
  const double t73 = P.t73;
  const double t74 = P.t74;
  const double t75 = P.t75;
  const double t76 = P.t76;
  const double t77 = P.t77;
  const double t78 = P.t78;
  return -13.0/2.0*m*t37 - m*t60 + (43.0/2.0)*m*t62 + 32*m*t70 + 18*m*t71 + (11.0/2.0)*m*t73 - t10*t5*v + t10*t6*t7 + t15*t16 - 4*t15 - t16*t57 + t17*t18 + (15.0/2.0)*t17*t21 - t19*t20 + t22*t27 - t22*t28 - t22 - 14*t24*v + (25.0/2.0)*t25*v - 32*t26*t57 - t28*t52 - t3*t5 + t3 + 10*t34*t35 + 4*t34*t41 - t36*t39 - 7.0/2.0*t36*t41 - t38*t8 + (5.0/2.0)*t4*t52 + (29.0/2.0)*t4*t57 - 8*t40*t41 + t42*t47 - t49*t58 + 12*t49 + t50*t51 - t50*t61 - t52*t58 + t52 - t53*t54 + t54*t55 + 6*t57 + (27.0/2.0)*t63*v + t64*t67 + 2*t64*t68 + 16*t66*t77 + (3.0/2.0)*t68*t75 + (17.0/2.0)*t68*t77 + 10*t72*v + (13.0/2.0)*t74*v + t75*t76 + t75*t78 + t76*t77;
}

inline double lognormal2_dexpected_var_var_var(const Lognormal2DexpectedPar& P) {
  const double t1 = P.t1;
  const double t11 = P.t11;
  const double t12 = P.t12;
  const double t34 = P.t34;
  const double t36 = P.t36;
  const double t40 = P.t40;
  const double t53 = P.t53;
  const double t61 = P.t61;
  const double t68 = P.t68;
  const double t86 = P.t86;
  const double t94 = P.t94;
  const double t95 = P.t95;
  const double t100 = P.t100;
  const double t101 = P.t101;
  const double t102 = P.t102;
  const double t103 = P.t103;
  const double t104 = P.t104;
  return t1*t104 + t1*t53 + t1*t94 + t100*t86 + t101*t11 + t101*t12 + t102*t11 + t102*t12 + t103*t11 + t103*t12 - t34 - 1.0/8.0*t36 - t40 + t61 + (1.0/4.0)*t68*t86 + t95;
}

inline void lognormal2_score_curv(int k, double y, const double* th,
                                  double* out) {
  const double m = th[0], v = th[1];
  const double C = std::log(y / m);
  const Lognormal2GradientPar Pg = lognormal2_gradient_par(m, v);
  const Lognormal2HessianPar Ph = lognormal2_hessian_par(m, v);
  if (k == 0) {
    out[0] = lognormal2_score_mean(Pg, C);
    out[1] = lognormal2_hess_mean_mean(Ph, C);
  } else {
    out[0] = lognormal2_score_var(Pg, C);
    out[1] = lognormal2_hess_var_var(Ph, C);
  }
}

inline void lognormal2_info_dinfo(int k, double y, const double* th,
                                  double* out) {
  const double m = th[0], v = th[1];
  const Lognormal2ExpectedPar Pe = lognormal2_expected_par(m, v);
  const Lognormal2DexpectedPar Pd = lognormal2_dexpected_par(m, v);
  if (k == 0) {
    out[0] = lognormal2_expected_mean_mean(Pe, v);
    out[1] = lognormal2_dexpected_mean_mean_mean(Pd, m, v);
  } else {
    out[0] = lognormal2_expected_var_var(Pe);
    out[1] = lognormal2_dexpected_var_var_var(Pd);
  }
}

} // namespace d7

#endif
