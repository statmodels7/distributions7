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
// The expected Hessian and its first and second derivatives (lognormal2_expected_all,
// lognormal2_dexpected1_all, lognormal2_dexpected2_all) take m and v and carry S and
// W = 1/(m^2 + v) as symbols; they have no y.

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

// E_ab = E[l_ab] on the parameter scale, components (mean_mean, var_var, mean_var).
// With S = log1p(v/m^2) and W = 1/(m^2 + v) carried as symbols, E = J' diag(-1/S, -1/(2 S^2)) J,
// J being the Jacobian of (mu, S) in (m, v). The expressions are the closed forms factored by
// sympy (misure_d34e/agents/lognormal2_low/gen.py), so that no power of m^2 + v is expanded.
inline void lognormal2_expected_all(double m, double v, double* o) {
  const double S = std::log1p(v / (m * m));
  const double W = 1.0 / (m * m + v);
  const double x0 = std::pow(S, -2);
  const double x1 = W*v;
  const double x2 = 2*x1;
  const double x3 = std::pow(W, 2);
  const double x4 = std::pow(v, 2)*x3;
  o[0] = -x0*(S*x2 + S*x4 + S + 2*x4)/std::pow(m, 2);
  o[1] = -1.0/4.0*x0*x3*(S + 2);
  o[2] = (1.0/2.0)*W*x0*(S*x1 + S + x2)/m;
}

// First derivatives of E_ab in the same S and W form, components (mean_mean_mean, mean_mean_var,
// var_var_mean, var_var_var, mean_var_mean, mean_var_var), the last index being the parameter
// differentiated.
inline void lognormal2_dexpected1_all(double m, double v, double* o) {
  const double S = std::log1p(v / (m * m));
  const double W = 1.0 / (m * m + v);
  const double x0 = std::pow(S, -3);
  const double x1 = std::pow(S, 2);
  const double x2 = 2*x1;
  const double x3 = W*v;
  const double x4 = x2*x3;
  const double x5 = std::pow(W, 3);
  const double x6 = std::pow(v, 3);
  const double x7 = std::pow(W, 2);
  const double x8 = std::pow(v, 2);
  const double x9 = x7*x8;
  const double x10 = std::pow(m, 2);
  const double x11 = 4*S*x10;
  const double x12 = x5*x8;
  const double x13 = x10*x2;
  const double x14 = S*v;
  const double x15 = W*x14;
  const double x16 = 5*S;
  const double x17 = W*x0/x10;
  const double x18 = x0*x5;
  const double x19 = (1.0/2.0)/m;
  const double x20 = x10*x7;
  const double x21 = 4*v;
  o[0] = -2*x0*(S*W*v + S*x5*x6 - v*x13*x7 - x1*x9 - x1 - x11*x12 - x12*x13 - x4 + 4*x5*x6)/std::pow(m, 3);
  o[1] = x17*(S - 2*x15 + x16*x9 + x2*x9 - x2 + 4*x9);
  o[2] = -x18*x19*(S*v + 4*v - x11 - x13);
  o[3] = (1.0/4.0)*x18*(x16 + x2 + 4);
  o[4] = (1.0/2.0)*x17*(2*S*x7*x8 - W*x13 - x1*x20*x21 - x1*x3 - x1 - 8*x14*x20 + 8*x7*x8);
  o[5] = -x0*x19*x7*(-S + W*x21 + 5*x15 + x4);
}

// Second derivatives of E_ab in the same form, components (mean_mean_mean_mean, mean_mean_var_var,
// mean_mean_mean_var, var_var_mean_mean, var_var_var_var, var_var_mean_var, mean_var_mean_mean,
// mean_var_var_var, mean_var_mean_var), the last two indices being the parameters differentiated.
inline void lognormal2_dexpected2_all(double m, double v, double* o) {
  const double S = std::log1p(v / (m * m));
  const double W = 1.0 / (m * m + v);
  const double x0 = std::pow(m, 4);
  const double x1 = std::pow(S, 3);
  const double x2 = 3*x1;
  const double x3 = std::pow(S, 2);
  const double x4 = v*x3;
  const double x5 = W*x4;
  const double x6 = 5*x5;
  const double x7 = 6*x1;
  const double x8 = W*v;
  const double x9 = x7*x8;
  const double x10 = std::pow(W, 4);
  const double x11 = std::pow(v, 4)*x10;
  const double x12 = 4*S;
  const double x13 = std::pow(W, 2);
  const double x14 = std::pow(v, 2);
  const double x15 = x13*x14;
  const double x16 = std::pow(W, 3);
  const double x17 = std::pow(v, 3);
  const double x18 = x16*x17;
  const double x19 = 4*x3;
  const double x20 = x15*x19;
  const double x21 = x10*x17;
  const double x22 = std::pow(m, 2);
  const double x23 = 40*x22;
  const double x24 = S*x23;
  const double x25 = 2*x3;
  const double x26 = x13*x22;
  const double x27 = v*x26;
  const double x28 = -x25*x27;
  const double x29 = x0*x1;
  const double x30 = v*x16;
  const double x31 = x22*x3;
  const double x32 = 10*x31;
  const double x33 = 24*x0;
  const double x34 = x10*x14;
  const double x35 = x14*x16;
  const double x36 = x22*x35;
  const double x37 = x36*x7;
  const double x38 = 12*x0*x1;
  const double x39 = std::pow(S, -4);
  const double x40 = 2*x39;
  const double x41 = 2*x1;
  const double x42 = -x41;
  const double x43 = 12*x8;
  const double x44 = x1*x8;
  const double x45 = 12*x15;
  const double x46 = 22*S;
  const double x47 = 17*x3;
  const double x48 = x39/x22;
  const double x49 = x13*x48;
  const double x50 = 8*S;
  const double x51 = W/std::pow(m, 3);
  const double x52 = x22*x41;
  const double x53 = S*x8;
  const double x54 = (1.0/2.0)*x16;
  const double x55 = W*x22;
  const double x56 = x39*x54/m;
  o[0] = -x40*(-12*S*x18 + x11*x12 + 24*x11 + x12*x15 + x15*x2 - 5*x18*x3 + x2 - x20 - x21*x24 - x21*x32 + x27*x7 + x28 + 8*x29*x30 + x3*x33*x34 + x34*x38 + x37 - x6 + x9)/x0;
  o[1] = -x49*(-S*x43 + 2*S + x15*x46 + x15*x47 + x15*x7 + x3 + x42 - 4*x44 + x45 - 14*x5);
  o[2] = x40*x51*(2*S*W*v + 14*S*x16*x17 - S*x45 + 2*W*x1*x22 + 6*v*x13*x22*x3 - x15*x41 + 3*x16*x17*x3 + 12*x16*x17 - x20 - x25 - 14*x31*x35 - x36*x50 - x37 - x42 - x5);
  o[3] = -x48*x54*(4*S*W*x14 + 12*W*x0*x1 + 24*W*x0*x3 + 24*W*x14 - v*x12 - x19*x22 - 10*x22*x5 - x23*x53 - x4 - x52);
  o[4] = -1.0/4.0*x10*x39*(x46 + x47 + x7 + 12);
  o[5] = x56*(14*S*W*v + 3*W*v*x3 + 12*W*v - 14*W*x22*x3 - x12 - x3 - x50*x55 - x55*x7);
  o[6] = x39*x51*(W*x1*x22 + v*x13*x52 + x1 + x12*x18 + 4*x13*x29 - 3*x15*x3 - x15*x50 + x16*x33*x4 + 24*x18 - x24*x35 + x28 + x30*x38 - x32*x35 + x44 - x5);
  o[7] = x56*(-6*S - 7*x3 + x42 + x43 + x46*x8 + 17*x5 + x9);
  o[8] = -1.0/2.0*x49*(28*S*x13*x14 - 16*S*x27 + 6*W*x22*x3 - 12*x1*x27 + 6*x13*x14*x3 + 24*x13*x14 - 28*x26*x4 - x3 - x41*x8 - 16*x53 - x6);
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
  double e[3], d[6];
  lognormal2_expected_all(m, v, e);
  lognormal2_dexpected1_all(m, v, d);
  if (k == 0) {
    out[0] = e[0];
    out[1] = d[0];
  } else {
    out[0] = e[1];
    out[1] = d[3];
  }
}

// Mixed components of the observed and the expected Hessian, and of the
// derivatives of the expected Hessian, each copied from the kernel of
// lognormal2.cpp that writes it (the line of the kernel in the comment).
// The parameter-only quantities are read from the struct, as above.

// lognormal2_hessian_cpp, o_mean_var
inline double lognormal2_hess_mean_var(const Lognormal2HessianPar& P, double m, double C) {
  const double t2 = P.t2;
  const double t3 = P.t3;
  const double t4 = P.t4;
  const double t6 = P.t6;
  const double t12 = P.t12;
  const double t14 = P.t14;
  const double t15 = P.t15;
  const double t16 = P.t16;
  const double t5 = std::pow(C, 2);
  const double t7 = t5*t6;
  return -C*t15 + m*t12*t2 + m*t12*t3 + (1.0/4.0)*m*t12 + 2*t14*t4*t5*t6 - t15 - t16*t2*t5 - 2*t16*t7;
}

// The routers of the joint registry, one per order: every component at one
// observation, on the parameter scale. The Hessian and the expected
// information fill a full 2 x 2 matrix and the derivative of the expected
// information a 2 x 2 x 2 array, column-major, out[a + 2 b + 4 c] being the
// derivative of E[l_ab] in parameter c. The parameters are (mean, var).
inline void lognormal2_joint_grad(double y, const double* th, double* out) {
  const double m = th[0], v = th[1];
  const double C = std::log(y / m);
  const Lognormal2GradientPar Pg = lognormal2_gradient_par(m, v);
  out[0] = lognormal2_score_mean(Pg, C);
  out[1] = lognormal2_score_var(Pg, C);
}

inline void lognormal2_joint_hess(double y, const double* th, double* out) {
  const double m = th[0], v = th[1];
  const double C = std::log(y / m);
  const Lognormal2HessianPar Ph = lognormal2_hessian_par(m, v);
  out[0] = lognormal2_hess_mean_mean(Ph, C);
  out[3] = lognormal2_hess_var_var(Ph, C);
  out[1] = out[2] = lognormal2_hess_mean_var(Ph, m, C);
}

inline void lognormal2_joint_info(double y, const double* th, double* out) {
  double e[3];
  lognormal2_expected_all(th[0], th[1], e);
  out[0] = e[0];
  out[3] = e[1];
  out[1] = out[2] = e[2];
}

inline void lognormal2_joint_dinfo(double y, const double* th, double* out) {
  double d[6];
  lognormal2_dexpected1_all(th[0], th[1], d);
  for (int j = 0; j < 8; ++j) out[j] = 0.0;
  out[0] = d[0];
  out[4] = d[1];
  out[3] = d[2];
  out[7] = d[3];
  out[1] = out[2] = d[4];
  out[5] = out[6] = d[5];
}

} // namespace d7

#endif
