#ifndef D7_PT_LOGNORMAL2_H
#define D7_PT_LOGNORMAL2_H

#include <Rcpp.h>
#include <cmath>
#include "pt_fdb2.h"
#include "pt_parts.h"

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

// The observed second derivatives, and the response derivatives that involve the
// second parameter derivatives, in the variables r = v/m^2 and u = 1/(1 + r).
// With l = F(S, C) - log y + const, F_S = (C^2 - S - S^2/4)/(2 S^2), F_C = -C/S - 1/2,
// the partials of S are monomials in r and u,
//   S_m = -2 r u/m,  S_v = u/m^2,  S_mm = 2 r (3 + r) u^2/m^2,
//   S_mv = -2 u^2/m^3,  S_vv = -u^2/m^4,
// so that no difference of nearly equal terms is formed when r is small (the
// expanded closed forms lose a factor 1/r in S_m and S_mm, since m W - 1/m is
// -v/(m (m^2 + v))). C = log(y/m) is taken by log1p((y - m)/m) near y = m.
struct Lognormal2Base { double m, S, Sm, Sv, Smm, Smv, Svv, Cm, Cmm; };

inline Lognormal2Base lognormal2_base(double m, double v) {
  const double im = 1.0 / m, im2 = im * im;
  const double r = v * im2, u = 1.0 / (1.0 + r);
  Lognormal2Base B;
  B.m = m;
  B.S = std::log1p(r);
  B.Sm = -2.0 * r * u * im;
  B.Sv = u * im2;
  B.Smm = 2.0 * r * (3.0 + r) * u * u * im2;
  B.Smv = -2.0 * u * u * im2 * im;
  B.Svv = -u * u * im2 * im2;
  B.Cm = -im;
  B.Cmm = im2;
  return B;
}

// C = log(y/m): log1p near y = m, the difference of logarithms elsewhere
// (as in lognormal2_obs.cpp).
inline double lognormal2_logratio(double y, double m) {
  return std::fabs(y - m) < 0.5 * m ? std::log1p((y - m) / m)
                                    : std::log(y) - std::log(m);
}

inline double lognormal2_hess_mean_mean(const Lognormal2Base& B, double C) {
  const double S = B.S, iS = 1.0 / S;
  const double FS = 0.5 * (C * C - S - 0.25 * S * S) * iS * iS;
  const double FSS = 0.5 * (S - 2.0 * C * C) * iS * iS * iS;
  const double FC = -C * iS - 0.5;
  return FSS * B.Sm * B.Sm + FS * B.Smm + 2.0 * (C * iS * iS) * B.Sm * B.Cm
         - iS * B.Cm * B.Cm + FC * B.Cmm;
}

inline double lognormal2_hess_mean_var(const Lognormal2Base& B, double C) {
  const double S = B.S, iS = 1.0 / S;
  const double FS = 0.5 * (C * C - S - 0.25 * S * S) * iS * iS;
  const double FSS = 0.5 * (S - 2.0 * C * C) * iS * iS * iS;
  return FSS * B.Sm * B.Sv + FS * B.Smv + (C * iS * iS) * B.Cm * B.Sv;
}

inline double lognormal2_hess_var_var(const Lognormal2Base& B, double C) {
  const double S = B.S, iS = 1.0 / S;
  const double FS = 0.5 * (C * C - S - 0.25 * S * S) * iS * iS;
  const double FSS = 0.5 * (S - 2.0 * C * C) * iS * iS * iS;
  return FSS * B.Sv * B.Sv + FS * B.Svv;
}

// d^2/dy d theta l = (1/y) d_theta F_C, and d^2/dy^2 d theta l = (d_theta F_CC - d_theta F_C)/y^2
inline double lognormal2_cross_y_mean(const Lognormal2Base& B, double y, double C) {
  const double iS = 1.0 / B.S;
  return (C * iS * iS * B.Sm - iS * B.Cm) / y;
}
inline double lognormal2_cross_y_var(const Lognormal2Base& B, double y, double C) {
  const double iS = 1.0 / B.S;
  return C * iS * iS * B.Sv / y;
}
inline double lognormal2_cross2_y_mean(const Lognormal2Base& B, double y, double C) {
  const double iS = 1.0 / B.S;
  return (iS * iS * B.Sm - (C * iS * iS * B.Sm - iS * B.Cm)) / (y * y);
}
inline double lognormal2_cross2_y_var(const Lognormal2Base& B, double y, double C) {
  const double iS = 1.0 / B.S;
  return (iS * iS * B.Sv - C * iS * iS * B.Sv) / (y * y);
}

// (1/y) d_a d_b F_C, with G = F_C: G_S = C/S^2, G_C = -1/S, G_SS = -2C/S^3, G_SC = 1/S^2.
inline void lognormal2_grad_y_hess_all(const Lognormal2Base& B, double y, double C, double* o) {
  const double iS = 1.0 / B.S, iS2 = iS * iS;
  const double GS = C * iS2, GC = -iS, GSS = -2.0 * C * iS2 * iS, GSC = iS2;
  o[0] = (GSS * B.Sm * B.Sm + GS * B.Smm + 2.0 * GSC * B.Sm * B.Cm + GC * B.Cmm) / y;
  o[1] = (GSS * B.Sm * B.Sv + GS * B.Smv + GSC * B.Cm * B.Sv) / y;
  o[2] = (GSS * B.Sv * B.Sv + GS * B.Svv) / y;
}

// (d_a d_b F_CC - d_a d_b F_C)/y^2, with F_CC = -1/S: H_S = 1/S^2, H_SS = -2/S^3.
inline void lognormal2_hess_y_hess_all(const Lognormal2Base& B, double y, double C, double* o) {
  const double iS = 1.0 / B.S, iS2 = iS * iS;
  const double GS = C * iS2, GC = -iS, GSS = -2.0 * C * iS2 * iS, GSC = iS2;
  const double HS = iS2, HSS = -2.0 * iS2 * iS;
  const double y2 = y * y;
  o[0] = (HSS * B.Sm * B.Sm + HS * B.Smm
          - (GSS * B.Sm * B.Sm + GS * B.Smm + 2.0 * GSC * B.Sm * B.Cm + GC * B.Cmm)) / y2;
  o[1] = (HSS * B.Sm * B.Sv + HS * B.Smv - (GSS * B.Sm * B.Sv + GS * B.Smv + GSC * B.Cm * B.Sv)) / y2;
  o[2] = (HSS * B.Sv * B.Sv + HS * B.Svv - (GSS * B.Sv * B.Sv + GS * B.Svv)) / y2;
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
  const Lognormal2Base B = lognormal2_base(m, v);
  const double C1 = lognormal2_logratio(y, m);
  if (k == 0) {
    out[0] = lognormal2_score_mean(Pg, C);
    out[1] = lognormal2_hess_mean_mean(B, C1);
  } else {
    out[0] = lognormal2_score_var(Pg, C);
    out[1] = lognormal2_hess_var_var(B, C1);
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
  const double C = lognormal2_logratio(y, m);
  const Lognormal2Base B = lognormal2_base(m, v);
  out[0] = lognormal2_hess_mean_mean(B, C);
  out[3] = lognormal2_hess_var_var(B, C);
  out[1] = out[2] = lognormal2_hess_mean_var(B, C);
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

// The parts route (pt_parts.h). The costly quantities are the logs: log(y/m)
// and the logs of 1 + v/m^2 (S). The Hessian and the gradient reach them
// through different expressions (joint_grad forms log(y/m) and log1p(v/(m m));
// joint_hess forms lognormal2_logratio and log1p(v (1/m)^2)), and each route
// keeps its own pair so that the orders agree bit for bit with their joints.
// Orders three and four take the S of ln2_par, log1p(v (1/m) (1/m)), and the
// parameter tables of Faa di Bruno that lognormal2_obs.cpp builds from it.
// Index tri(i, j) = (i + j)(i + j + 1)/2 + j, with i the mean and j the
// variance order, i + j <= 4 (15 entries):
//   [0] level                        always
//   [1] y                            observed parts; NaN in the expected parts
//   [2] m                            always
//   [3] v                            always
//   [4] C_g = log(y/m)               level >= 1 (the gradient's C)
//   [5] S_g = log1p(v/(m m))         level >= 1 (the gradient's S, the expected information's S)
//   [6] C_h = lognormal2_logratio    level >= 2 (the Hessian's C, orders 3 and 4)
//   [7] S_b = log1p(v (1/m)^2)       level >= 2 (the Hessian's S, lognormal2_base)
//   [8] S_k = log1p(v (1/m) (1/m))   level >= 3 (the S of ln2_par)
//   [9 + tri(i, j)]  dH2[i][j], i + j <= level: level 3 for i + j <= 3, level 4 for i + j = 4
//   [24 + tri(i, j)] dH1[i][j], 3 <= i + j <= level: level 3 for i + j = 3, level 4 for i + j = 4
// Entries above the level are NaN. The expected parts keep [5] only.
inline int lognormal2_parts_tri(int i, int j) {
    const int t = i + j;
    return t * (t + 1) / 2 + j;
}

// The gradient's parameter quantities of lognormal2_gradient_par, with S given.
inline Lognormal2GradientPar lognormal2_gradient_par_from(double m, double v, double S) {
    const double t0 = 1.0 / m;
    const double t1 = 1.0 / S;
    const double t2 = t0 * t1;
    const double t3 = 1.0 / (std::pow(m, 2) + v);
    const double t4 = m * t3;
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

// The Hessian's quantities of lognormal2_base, with S given.
inline Lognormal2Base lognormal2_base_from(double m, double v, double S) {
    const double im = 1.0 / m, im2 = im * im;
    const double r = v * im2, u = 1.0 / (1.0 + r);
    Lognormal2Base B;
    B.m = m;
    B.S = S;
    B.Sm = -2.0 * r * u * im;
    B.Sv = u * im2;
    B.Smm = 2.0 * r * (3.0 + r) * u * u * im2;
    B.Smv = -2.0 * u * u * im2 * im;
    B.Svv = -u * u * im2 * im2;
    B.Cm = -im;
    B.Cmm = im2;
    return B;
}

// lognormal2_expected_all with S given.
inline void lognormal2_expected_from(double m, double v, double S, double* o) {
    const double W = 1.0 / (m * m + v);
    const double x0 = std::pow(S, -2);
    const double x1 = W * v;
    const double x2 = 2 * x1;
    const double x3 = std::pow(W, 2);
    const double x4 = std::pow(v, 2) * x3;
    o[0] = -x0 * (S * x2 + S * x4 + S + 2 * x4) / std::pow(m, 2);
    o[1] = -1.0 / 4.0 * x0 * x3 * (S + 2);
    o[2] = (1.0 / 2.0) * W * x0 * (S * x1 + S + x2) / m;
}

// C^(k) = (-1)^k (k-1)! / m^k for k = 1..K, as ln2_par forms it (Cm[0] unused).
inline void lognormal2_parts_cm(double m, int K, double* Cm) {
    const double im = 1.0 / m;
    double c = -im;
    Cm[0] = 0.0;
    for (int k = 1; k <= K; ++k) { Cm[k] = c; c *= -k * im; }
}

// (C^2)^(k) by Leibniz, k = 0..K, as ln2_c2 forms it.
inline void lognormal2_parts_c2(double C, const double* Cm, int K, double* C2d) {
    double Cd[6];
    Cd[0] = C;
    for (int k = 1; k <= K; ++k) Cd[k] = Cm[k];
    for (int k = 0; k <= K; ++k) {
        double s = 0.0, binom = 1.0;
        for (int l = 0; l <= k; ++l) {
            s += binom * Cd[l] * Cd[k - l];
            binom = binom * (k - l) / (l + 1);
        }
        C2d[k] = s;
    }
}

// The Faa di Bruno tables dH2 (all i + j <= K) and dH1 (3 <= i + j <= K) of
// ln2_par, written into the buffer. The tables of order K agree with those of
// the kernel at any K >= i + j, because each entry depends only on its own
// total order.
inline void lognormal2_parts_tables(double m, double v, int K, double* b) {
    const double im = 1.0 / m, r = v * im * im;
    double dr[6][6] = {{0.0}};
    double fm = im * im;
    for (int p = 0; p <= K; ++p) {
        dr[p][0] = v * fm;
        dr[p][1] = fm;
        fm *= -(p + 2.0) * im;
    }
    double sd[6] = {0.0};
    {
        const double i1 = 1.0 / (1.0 + r);
        double f = i1;
        for (int k = 1; k <= K; ++k) { sd[k] = f; f *= -k * i1; }
    }
    const double S = b[8];
    double dS[6][6] = {{0.0}};
    dS[0][0] = S;
    for (int p = 0; p <= K; ++p)
        for (int q = 0; p + q <= K; ++q)
            if (p + q >= 1) dS[p][q] = d7::fdb2_compose(p, q, sd, dr);
    double h1d[6] = {0.0}, h2d[6] = {0.0};
    {
        const double iS = 1.0 / S;
        double f = -0.5 * iS;
        double g = 0.5 * iS * iS;
        for (int k = 1; k <= K; ++k) {
            h1d[k] = f;
            h2d[k] = g;
            f *= -k * iS;
            g *= -(k + 1.0) * iS;
        }
        h1d[1] -= 0.125;
    }
    for (int p = 0; p <= K; ++p)
        for (int q = 0; p + q <= K; ++q) {
            if (p + q == 0) {
                b[9] = -0.5 / S;
                continue;
            }
            b[9 + lognormal2_parts_tri(p, q)] = d7::fdb2_compose(p, q, h2d, dS);
            if (p + q >= 3)
                b[24 + lognormal2_parts_tri(p, q)] = d7::fdb2_compose(p, q, h1d, dS);
        }
}

// The (i, j) component of ln2_component, read from the buffer.
inline double lognormal2_parts_component(int i, int j, const double* b,
                                         const double* Cm, const double* C2d) {
    double s = b[24 + lognormal2_parts_tri(i, j)];
    double binom = 1.0;
    for (int kk = 0; kk <= i; ++kk) {
        s += binom * b[9 + lognormal2_parts_tri(i - kk, j)] * C2d[kk];
        binom = binom * (i - kk) / (kk + 1);
    }
    if (j == 0) s -= 0.5 * Cm[i];
    return s;
}

inline void lognormal2_parts(double y, const double* th, int level, double* b) {
    const double m = th[0], v = th[1];
    b[0] = level; b[1] = y; b[2] = m; b[3] = v;
    for (int j = 4; j < 39; ++j) b[j] = R_NaN;
    if (level >= 1) {
        b[4] = std::log(y / m);
        b[5] = std::log1p(v / (m * m));
    }
    if (level >= 2) {
        const double im = 1.0 / m, im2 = im * im;
        b[6] = lognormal2_logratio(y, m);
        b[7] = std::log1p(v * im2);
    }
    if (level >= 3) {
        const double im = 1.0 / m;
        b[8] = std::log1p(v * im * im);
        lognormal2_parts_tables(m, v, level, b);
    }
}

inline void lognormal2_parts_expected(const double* th, int level, double* b) {
    b[0] = level; b[1] = R_NaN; b[2] = th[0]; b[3] = th[1]; b[4] = R_NaN;
    b[5] = level >= 2 ? std::log1p(th[1] / (th[0] * th[0])) : R_NaN;
}

inline void lognormal2_parts_grad(const double* b, double* out) {
    const Lognormal2GradientPar Pg = lognormal2_gradient_par_from(b[2], b[3], b[5]);
    out[0] = lognormal2_score_mean(Pg, b[4]);
    out[1] = lognormal2_score_var(Pg, b[4]);
}

inline void lognormal2_parts_hess(const double* b, double* out) {
    const Lognormal2Base B = lognormal2_base_from(b[2], b[3], b[7]);
    out[0] = lognormal2_hess_mean_mean(B, b[6]);
    out[3] = lognormal2_hess_var_var(B, b[6]);
    out[1] = out[2] = lognormal2_hess_mean_var(B, b[6]);
}

// from the expected parts
inline void lognormal2_parts_info(const double* b, double* out) {
    double e[3];
    lognormal2_expected_from(b[2], b[3], b[5], e);
    out[0] = e[0];
    out[3] = e[1];
    out[1] = out[2] = e[2];
}

// third order: lognormal2_deriv3_cpp (lognormal2_obs.cpp), unique components in
// the sorted order mean mean mean, mean mean var, mean var var, var var var
inline void lognormal2_parts_d3(const double* b, double* out) {
    double Cm[6], C2d[6];
    lognormal2_parts_cm(b[2], 3, Cm);
    lognormal2_parts_c2(b[6], Cm, 3, C2d);
    double u[4];
    u[0] = lognormal2_parts_component(3, 0, b, Cm, C2d);
    u[1] = lognormal2_parts_component(2, 1, b, Cm, C2d);
    u[2] = lognormal2_parts_component(1, 2, b, Cm, C2d);
    u[3] = lognormal2_parts_component(0, 3, b, Cm, C2d);
    sym_fill(2, 3, u, out);
}

// fourth order: lognormal2_deriv4_cpp (lognormal2_obs.cpp), unique components in
// the sorted order mean^4, mean^3 var, mean^2 var^2, mean var^3, var^4
inline void lognormal2_parts_d4(const double* b, double* out) {
    double Cm[6], C2d[6];
    lognormal2_parts_cm(b[2], 4, Cm);
    lognormal2_parts_c2(b[6], Cm, 4, C2d);
    double u[5];
    u[0] = lognormal2_parts_component(4, 0, b, Cm, C2d);
    u[1] = lognormal2_parts_component(3, 1, b, Cm, C2d);
    u[2] = lognormal2_parts_component(2, 2, b, Cm, C2d);
    u[3] = lognormal2_parts_component(1, 3, b, Cm, C2d);
    u[4] = lognormal2_parts_component(0, 4, b, Cm, C2d);
    sym_fill(2, 4, u, out);
}

} // namespace d7

#endif
