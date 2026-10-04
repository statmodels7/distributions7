#ifndef D7_PT_SKEWNORMAL2_H
#define D7_PT_SKEWNORMAL2_H

#include <Rcpp.h>
#include <cmath>
#include "pt_loc_scale.h"
#include "pt_sqrt.h"

// The skew normal in its centered parametrization (mu, sigma, gamma1); see
// skewnormal2.cpp for the construction. With w = (y - mu)/sigma,
// r = sign(g) (|g|/c)^(1/3), Dq = sqrt(b^2 - (1 - b^2) r^2),
// X = r (w + r) / (Dq sqrt(1 + r^2)) and Z1 = phi(X)/Phi(X), every component
// is a combination of G_ak = d^a_w d^k_g F. One function per G_ak and per
// component for the score, the diagonal of the Hessian and the diagonal of
// the third derivatives, called by the vector kernels (skewnormal2.cpp), by
// the quadrature of the expected information and by the scalar registry.
// The G_ak are the generated kernels' own expressions, lifted with their
// temporaries; G_10 comes in two forms because the score's kernel and the
// Hessian's write its series branch differently.

#define SN2_XC 0.4
#define SN2_RC 0.4
// |g| below which the expected information is a series; sn2_ge() in R
#define SN2_GE 3e-3

namespace d7 {

const double SN2_C = (4.0 - M_PI) / 2.0;
const double SN2_B2 = 2.0 / M_PI;

inline double sn2_r(double g) {
  if (g == 0.0) return 0.0;
  const double a = std::cbrt(std::fabs(g) / SN2_C);
  return g > 0.0 ? a : -a;
}

// sqrt(b^2 - (1 - b^2) r^2) = b sqrt(1 - delta^2), with
// 1 - delta^2 = -expm1((2/3) log(|g|/g_max)): the difference reaches 1e-16
// at the ceiling of the skewness, where this keeps its digits
inline double sn2_Dq(double g) {
  const double gmax = (4.0 - M_PI) / 2.0 *
    std::pow(d7::sqrt_cr(SN2_B2) / d7::sqrt_cr(1.0 - SN2_B2), 3.0);
  if (g == 0.0) return d7::sqrt_cr(SN2_B2);
  return d7::sqrt_cr(SN2_B2 * -std::expm1((2.0 / 3.0) * std::log(std::fabs(g) / gmax)));
}

// phi(x)/Phi(x) on the log scale, finite far below where both underflow
inline double sn2_zeta1(double x) {
  return std::exp(R::dnorm4(x, 0.0, 1.0, 1) - R::pnorm5(x, 0.0, 1.0, 1, 1));
}

// the expected information near g = 0: E_ij = s^-p sum_n e_n r^n, n = 0..36
// (skewnormal2.cpp)
extern const double SN2_E[6][37];
extern const int SN2_EP[6];

// sum_n e_n prod_(i<k) (n - 3i)/(3c) r^(n - 3k): the k-th derivative in g
inline double sn2_eser(const double* e, double r, int k) {
  double out = 0.0;
  for (int nn = 0; nn <= 36; nn++) {
    if (e[nn] == 0.0) continue;
    double fac = 1.0;
    for (int i = 0; i < k; i++) fac *= (nn - 3.0 * i) / (3.0 * SN2_C);
    if (fac == 0.0) continue;
    out += e[nn] * fac * std::pow(r, nn - 3 * k);
  }
  return out;
}

// d^j/ds^j s^-p = (-1)^j p (p+1) ... (p+j-1) s^(-p-j)
inline double sn2_spow(double s, int p, int j) {
  double c = 1.0;
  for (int i = 0; i < j; i++) c *= -(p + i);
  return c * std::pow(s, -p - j);
}

// G_ak from the series in r, k >= 1 (skewnormal2.cpp)
double sn2_G_series(int a, int k, double w, double r);

// where the kernels take G_ak, k >= 1, from the series rather than from the
// closed form
inline bool sn2_series_region(double X, double r) {
  return std::fabs(X) < SN2_XC && std::fabs(r) < SN2_RC;
}

// the shape-only temporaries of the G_ak below, one field per distinct
// expression, and one builder per set of G_ak a caller reads; a field is
// computed in the order the generated kernels define it
struct Sn2ShapeGrad {
  double r, Dq;
  double p0, p1, p2, p3, p4, p5;
};

inline Sn2ShapeGrad sn2_shape_grad(double g) {
  Sn2ShapeGrad P;
  P.r = sn2_r(g);
  P.Dq = sn2_Dq(g);
  P.p0 = std::pow(P.r, 2) + 1;
  P.p1 = std::pow(P.r, 2);
  P.p2 = P.p1 + 1;
  P.p3 = 1.0/P.p2;
  P.p4 = 1.0/P.Dq;
  P.p5 = std::pow(P.p2, -1.0/2.0);
  return P;
}

struct Sn2ShapeHess {
  double r, Dq;
  double p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17;
};

inline Sn2ShapeHess sn2_shape_hess(double g) {
  Sn2ShapeHess P;
  P.r = sn2_r(g);
  P.Dq = sn2_Dq(g);
  P.p1 = std::pow(P.r, 2);
  P.p2 = P.p1 + 1;
  P.p3 = 1.0/P.p2;
  P.p4 = 1.0/P.Dq;
  P.p5 = std::pow(P.p2, -1.0/2.0);
  P.p6 = std::pow(P.Dq, -2);
  P.p7 = P.p1*P.p6;
  P.p8 = 1.0/P.r;
  P.p9 = 1.0/(4 - M_PI);
  P.p10 = 1 - 2/M_PI;
  P.p11 = 2*P.r;
  P.p12 = std::pow(P.p2, -2);
  P.p13 = P.p1*P.p3;
  P.p14 = -P.p10;
  P.p15 = P.p14*P.p7;
  P.p16 = std::pow(P.r, 3);
  P.p17 = P.p14*P.p6;
  return P;
}

struct Sn2ShapeD3 {
  double r, Dq;
  double p1, p2, p3, p4, p5, p6, p7, p8, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37;
};

inline Sn2ShapeD3 sn2_shape_d3(double g) {
  Sn2ShapeD3 P;
  P.r = sn2_r(g);
  P.Dq = sn2_Dq(g);
  P.p1 = std::pow(P.r, 2);
  P.p2 = P.p1 + 1;
  P.p3 = 1.0/P.p2;
  P.p4 = 1.0/P.Dq;
  P.p5 = std::pow(P.p2, -1.0/2.0);
  P.p6 = std::pow(P.Dq, -2);
  P.p7 = P.p1*P.p6;
  P.p8 = 1.0/P.r;
  P.p10 = 1 - 2/M_PI;
  P.p11 = 2*P.r;
  P.p12 = std::pow(P.p2, -2);
  P.p13 = P.p1*P.p3;
  P.p14 = -P.p10;
  P.p15 = P.p14*P.p7;
  P.p16 = std::pow(P.r, 3);
  P.p17 = P.p14*P.p6;
  P.p18 = P.p10*P.p7;
  P.p19 = P.p13 + P.p15 - 1;
  P.p20 = 4 - M_PI;
  P.p21 = (2.0/3.0)/P.p20;
  P.p22 = std::pow(P.Dq, -3);
  P.p23 = 2*P.p12;
  P.p24 = -P.p1*P.p3 + 1;
  P.p25 = 2*P.p3;
  P.p26 = 1.0/P.p1;
  P.p27 = std::pow(P.p20, -2);
  P.p28 = 2*P.p26;
  P.p29 = P.p10*P.p6;
  P.p30 = 3*P.p12;
  P.p31 = P.p16*P.p17;
  P.p32 = 3*P.p18;
  P.p33 = std::pow(P.p2, -3);
  P.p34 = 8*P.p33;
  P.p35 = 2*P.p29;
  P.p36 = std::pow(P.r, 4);
  P.p37 = P.p3*P.r;
  return P;
}

struct Sn2ShapeAll {
  double r, Dq;
  double p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p20, p22, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37;
};

inline Sn2ShapeAll sn2_shape_all(double g) {
  Sn2ShapeAll P;
  P.r = sn2_r(g);
  P.Dq = sn2_Dq(g);
  P.p0 = std::pow(P.r, 2) + 1;
  P.p1 = std::pow(P.r, 2);
  P.p2 = P.p1 + 1;
  P.p3 = 1.0/P.p2;
  P.p4 = 1.0/P.Dq;
  P.p5 = std::pow(P.p2, -1.0/2.0);
  P.p6 = std::pow(P.Dq, -2);
  P.p7 = P.p1*P.p6;
  P.p8 = 1.0/P.r;
  P.p9 = 1.0/(4 - M_PI);
  P.p10 = 1 - 2/M_PI;
  P.p11 = 2*P.r;
  P.p12 = std::pow(P.p2, -2);
  P.p13 = P.p1*P.p3;
  P.p14 = -P.p10;
  P.p15 = P.p14*P.p7;
  P.p16 = std::pow(P.r, 3);
  P.p17 = P.p14*P.p6;
  P.p18 = P.p10*P.p7;
  P.p20 = 4 - M_PI;
  P.p22 = std::pow(P.Dq, -3);
  P.p25 = 2*P.p3;
  P.p26 = 1.0/P.p1;
  P.p27 = std::pow(P.p20, -2);
  P.p28 = 2*P.p26;
  P.p29 = P.p10*P.p6;
  P.p30 = 3*P.p12;
  P.p31 = P.p16*P.p17;
  P.p32 = 3*P.p18;
  P.p33 = std::pow(P.p2, -3);
  P.p34 = 8*P.p33;
  P.p35 = 2*P.p29;
  P.p36 = std::pow(P.r, 4);
  P.p37 = P.p3*P.r;
  return P;
}

template <class P_t>
inline double sn2_grad_G1_0(double w, double X, double Z1, const P_t& P) {
  (void) w; (void) X; (void) Z1;
  if (sn2_series_region(X, P.r)) {
    return -(P.r + w)/P.p0 + Z1*P.r/(P.Dq*d7::sqrt_cr(P.p0));
  }
  const double f0 = P.r + w;
  const double f4 = f0*P.p3;
  return Z1*P.p4*P.p5*P.r - f4;
}

template <class P_t>
inline double sn2_G1_0(double w, double X, double Z1, const P_t& P) {
  (void) w; (void) X; (void) Z1;
  if (sn2_series_region(X, P.r)) {
    return -P.p3*(P.r + w) + Z1*P.r/(P.Dq*d7::sqrt_cr(P.p2));
  }
  const double f0 = P.r + w;
  const double f4 = f0*P.p3;
  return Z1*P.p4*P.p5*P.r - f4;
}

template <class P_t>
inline double sn2_G2_0(double w, double X, double Z1, const P_t& P) {
  (void) w; (void) X; (void) Z1;
  if (sn2_series_region(X, P.r)) {
    return -P.p3*(1 + Z1*P.p1*(X + Z1)/std::pow(P.Dq, 2));
  }
  const double f9 = X + Z1;
  const double f10 = Z1*f9;
  return -P.p3*(f10*P.p7 + 1);
}

template <class P_t>
inline double sn2_G3_0(double w, double X, double Z1, const P_t& P) {
  (void) w; (void) X; (void) Z1;
  if (sn2_series_region(X, P.r)) {
    return Z1*std::pow(P.r, 3)*(std::pow(X, 2) + 3*X*Z1 + 2*std::pow(Z1, 2) - 1)/(std::pow(P.Dq, 3)*std::pow(P.p2, 3.0/2.0));
  }
  const double f13 = Z1*P.p22;
  const double f14 = f13*(std::pow(X, 2) + 3*X*Z1 + 2*std::pow(Z1, 2) - 1)/std::pow(P.p2, 3.0/2.0);
  return P.p16*f14;
}

// gs: G_ak from the series, read where sn2_series_region() holds
template <class P_t>
inline double sn2_G0_1(double w, double X, double Z1, const P_t& P,
                        double gs) {
  (void) w; (void) X; (void) Z1;
  if (sn2_series_region(X, P.r)) return gs;
  const double f0 = P.r + w;
  const double f4 = f0*P.p3;
  return -1.0/3.0*(Z1*P.p4*P.p5*(P.p1*f4 - 2*P.r - w + f0*P.p1*(-1 + M_2_PI)/std::pow(P.Dq, 2)) - std::pow(f0, 2)*P.r/std::pow(P.p2, 2) + P.p3*P.r + f4)/(P.p1*(2 - 1.0/2.0*M_PI));
}

// gs: G_ak from the series, read where sn2_series_region() holds
template <class P_t>
inline double sn2_G0_2(double w, double X, double Z1, const P_t& P,
                        double gs) {
  (void) w; (void) X; (void) Z1;
  if (sn2_series_region(X, P.r)) return gs;
  const double f0 = P.r + w;
  const double f4 = f0*P.p3;
  const double f9 = X + Z1;
  const double f10 = Z1*f9;
  const double f14 = f0*P.p10*P.p7;
  const double f16 = -f0*P.p1*P.p3 + P.p11 + w;
  const double f18 = f0*P.p12;
  const double f22 = Z1*P.p5;
  const double f23 = f22*P.p4;
  const double f25 = f0*P.p15 + P.p1*f4 - P.p11 - w;
  const double f26 = std::pow(f0, 2);
  return -2.0/9.0*P.p9*(f10*P.p8*std::pow(f25, 2)*P.p3*P.p6 - P.p8*(4*f0*P.p12*P.r + 2*P.p1*P.p12 - 4*P.p1*f26/std::pow(P.p2, 3) + P.p12*f26 - f23*(f0*P.p11*P.p17 - 3*f18*P.p16 + 2*P.p13 + P.p15 - P.p16*P.p17*f4 + 3*f4*P.r - 2) - 2*P.p3) + 2*(P.p12*f26*P.r - f23*f25 - P.p3*P.r - f4)/P.p1 + P.p10*f22*(-3*f14 - f16)/std::pow(P.Dq, 3))/(P.p16*(2 - 1.0/2.0*M_PI));
}

// gs: G_ak from the series, read where sn2_series_region() holds
template <class P_t>
inline double sn2_G0_3(double w, double X, double Z1, const P_t& P,
                        double gs) {
  (void) w; (void) X; (void) Z1;
  if (sn2_series_region(X, P.r)) return gs;
  const double f0 = P.r + w;
  const double f4 = f0*P.p3;
  const double f9 = X + Z1;
  const double f10 = Z1*f9;
  const double f13 = Z1*P.p22;
  const double f14 = f13*(std::pow(X, 2) + 3*X*Z1 + 2*std::pow(Z1, 2) - 1)/std::pow(P.p2, 3.0/2.0);
  const double f18 = f0*P.p18;
  const double f20 = -f0*P.p1*P.p3 + P.p11 + w;
  const double f23 = f0*P.p12;
  const double f28 = Z1*P.p5;
  const double f29 = f28*P.p4;
  const double f36 = f10*P.p6;
  const double f37 = P.p25*f36;
  const double f41 = f0*P.p15;
  const double f42 = P.p1*f4 - P.p11 - w;
  const double f43 = f41 + f42;
  const double f45 = -3*f18 - f20;
  const double f47 = f45*P.p29;
  const double f48 = 3*f4;
  const double f50 = f0*P.p30;
  const double f51 = P.p16*f50;
  const double f53 = f0*P.p17;
  const double f55 = P.p11*f53 + 2*P.p13 + P.p15 - f4*P.p31 + f48*P.r - f51 - 2;
  const double f56 = P.p8*f55;
  const double f57 = -f43*P.p28 + f47 + f56;
  const double f58 = P.p3*f36;
  const double f60 = P.p22*P.p10*f28;
  const double f63 = f0*P.p1;
  const double f66 = P.p8*f58;
  const double f67 = P.p5*(3*f0*P.p16*P.p10*P.p3*P.p6 + 3*f0*P.p3*P.r - 6*f0*P.p29*P.r + 2*P.p1*P.p3 - f51 - P.p32 - 2);
  const double f68 = f13*P.p10;
  const double f69 = f45*f60;
  const double f70 = std::pow(f43, 2);
  const double f72 = std::pow(f0, 2);
  const double f73 = P.p12*f72*P.r - f29*f43 - f4 - P.p37;
  const double f74 = P.p33*f72;
  const double f75 = 4*f0*P.p12*P.r + 2*P.p1*P.p12 - 4*P.p1*f74 + P.p12*f72 - f29*f55 - P.p25;
  return -4.0/27.0*P.p27*(f14*P.p26*std::pow(f43, 3) - P.p8*(P.p8*f69 + P.p8*(6*f0*P.p12 - P.p16*P.p34 + 24*P.p16*f72/std::pow(P.p2, 4) + 12*P.p12*P.r - f29*(15*f0*P.p33*P.p36 - 18*P.p1*f23 - 9*P.p16*P.p12 - 5*P.p15*f4 - P.p25*P.p31 + f48 + f50*P.p17*P.p36 + 4*P.p17*P.r + 2*f53 + 9*P.p37) - 24*P.p33*f63 - 12*f74*P.r) + f37*P.p26*f70 - f37*f43*f56 - 4*P.p26*f75 - f67*f68 + 6*f73/P.p16) + f43*f57*f66 + P.p28*(P.p8*f75 - P.p28*f73 - f66*f70 - f69) + f68*(P.p8*P.p25*f43*f45*P.p4*f9 + P.p8*f67 - P.p28*f45*P.p5 + f47*P.p5 + P.p5*P.p35*(6*f41 + f42)))/(P.p36*(2 - 1.0/2.0*M_PI));
}

// gs: G_ak from the series, read where sn2_series_region() holds
template <class P_t>
inline double sn2_hess_G1_1(double w, double X, double Z1, const P_t& P,
                             double gs) {
  (void) w; (void) X; (void) Z1;
  if (sn2_series_region(X, P.r)) return gs;
  const double f0 = P.r + w;
  const double f9 = X + Z1;
  const double f14 = f0*P.p10*P.p7;
  const double f16 = -f0*P.p1*P.p3 + P.p11 + w;
  const double f18 = f0*P.p12;
  const double f22 = Z1*P.p5;
  const double f23 = f22*P.p4;
  return (2.0/3.0)*P.p8*P.p9*(Z1*P.p3*P.p6*f9*(-f14 - f16) - P.p8*(-P.p11*f18 + f23*(P.p13 + P.p15 - 1) + P.p3));
}

// gs: G_ak from the series, read where sn2_series_region() holds
template <class P_t>
inline double sn2_G1_1(double w, double X, double Z1, const P_t& P,
                        double gs) {
  (void) w; (void) X; (void) Z1;
  if (sn2_series_region(X, P.r)) return gs;
  const double f0 = P.r + w;
  const double f9 = X + Z1;
  const double f18 = f0*P.p18;
  const double f20 = -f0*P.p1*P.p3 + P.p11 + w;
  const double f21 = -f18 - f20;
  const double f23 = f0*P.p12;
  const double f28 = Z1*P.p5;
  const double f29 = f28*P.p4;
  const double f30 = -P.p11*f23 + P.p19*f29 + P.p3;
  return P.p8*P.p21*(Z1*f21*P.p3*P.p6*f9 - P.p8*f30);
}

// gs: G_ak from the series, read where sn2_series_region() holds
template <class P_t>
inline double sn2_G2_1(double w, double X, double Z1, const P_t& P,
                        double gs) {
  (void) w; (void) X; (void) Z1;
  if (sn2_series_region(X, P.r)) return gs;
  const double f0 = P.r + w;
  const double f9 = X + Z1;
  const double f10 = Z1*f9;
  const double f13 = Z1*P.p22;
  const double f14 = f13*(std::pow(X, 2) + 3*X*Z1 + 2*std::pow(Z1, 2) - 1)/std::pow(P.p2, 3.0/2.0);
  const double f18 = f0*P.p18;
  const double f20 = -f0*P.p1*P.p3 + P.p11 + w;
  const double f21 = -f18 - f20;
  const double f36 = f10*P.p6;
  const double f37 = P.p25*f36;
  const double f38 = P.p8*f37;
  return P.p21*(-f14*f21 + P.p8*P.p23 + f38*(-P.p18 - P.p24));
}

// gs: G_ak from the series, read where sn2_series_region() holds
template <class P_t>
inline double sn2_G1_2(double w, double X, double Z1, const P_t& P,
                        double gs) {
  (void) w; (void) X; (void) Z1;
  if (sn2_series_region(X, P.r)) return gs;
  const double f0 = P.r + w;
  const double f4 = f0*P.p3;
  const double f9 = X + Z1;
  const double f10 = Z1*f9;
  const double f13 = Z1*P.p22;
  const double f14 = f13*(std::pow(X, 2) + 3*X*Z1 + 2*std::pow(Z1, 2) - 1)/std::pow(P.p2, 3.0/2.0);
  const double f18 = f0*P.p18;
  const double f20 = -f0*P.p1*P.p3 + P.p11 + w;
  const double f21 = -f18 - f20;
  const double f23 = f0*P.p12;
  const double f28 = Z1*P.p5;
  const double f29 = f28*P.p4;
  const double f30 = -P.p11*f23 + P.p19*f29 + P.p3;
  const double f36 = f10*P.p6;
  const double f37 = P.p25*f36;
  const double f38 = P.p8*f37;
  const double f41 = f0*P.p15;
  const double f42 = P.p1*f4 - P.p11 - w;
  const double f43 = f41 + f42;
  const double f45 = -3*f18 - f20;
  const double f47 = f45*P.p29;
  const double f48 = 3*f4;
  const double f50 = f0*P.p30;
  const double f51 = P.p16*f50;
  const double f53 = f0*P.p17;
  const double f55 = P.p11*f53 + 2*P.p13 + P.p15 - f4*P.p31 + f48*P.r - f51 - 2;
  const double f56 = P.p8*f55;
  const double f57 = -f43*P.p28 + f47 + f56;
  const double f58 = P.p3*f36;
  const double f60 = P.p22*P.p10*f28;
  const double f63 = f0*P.p1;
  return (4.0/9.0)*P.p26*P.p27*(f14*P.p8*std::pow(f21, 2) - P.p8*(-P.p8*(f0*P.p23 + 4*P.p12*P.r + f29*P.r*(P.p1*P.p30 - P.p13*P.p29 - 3*P.p3 + P.p35) - P.p34*f63) + P.p19*f38*f43 - f30*P.p28 + f60*(-P.p24 - P.p32)) + f57*f58);
}


// the components, from the G_ak
inline double skewnormal2_score_mu(double G1_0, double s) {
  const double u0 = 1.0/s;
  return -G1_0*u0;
}

inline double skewnormal2_score_sigma(double G1_0, double w, double s) {
  const double u0 = 1.0/s;
  return -u0*(G1_0*w + 1);
}

inline double skewnormal2_score_gamma1(double G0_1) {
  return G0_1;
}

inline double skewnormal2_hess_mu_mu(double G2_0, double s) {
  const double u0 = std::pow(s, -2);
  return G2_0*u0;
}

inline double skewnormal2_hess_sigma_sigma(double G1_0, double G2_0,
                                           double w, double s) {
  const double u0 = std::pow(s, -2);
  return u0*(2*G1_0*w + G2_0*std::pow(w, 2) + 1);
}

inline double skewnormal2_hess_gamma1_gamma1(double G0_2) {
  return G0_2;
}

inline double skewnormal2_d3_mu_mu_mu(double G3_0, double s) {
  const double u0 = std::pow(s, -3);
  return -G3_0*u0;
}

inline double skewnormal2_d3_sigma_sigma_sigma(double G1_0, double G2_0,
                                               double G3_0, double w,
                                               double s) {
  const double u0 = std::pow(s, -3);
  const double u2 = std::pow(w, 2);
  return -u0*(6*G1_0*w + 6*G2_0*u2 + G3_0*std::pow(w, 3) + 2);
}

inline double skewnormal2_d3_gamma1_gamma1_gamma1(double G0_3) {
  return G0_3;
}

// log f: skewnormal1's density at the direct parameters, mapped as
// sn_cp_to_dp() (skewnormal2_distrib.R) maps them, expression for
// expression, so that it is distrib_pdf()'s to the last bit: with
// c = sign(g) (2|g|/(4 - pi))^(1/3), xi = mu - s c,
// omega = s sqrt(1 + c^2), alpha = c / (b sqrt(1 - delta^2))
inline double skewnormal2_logpdf(double y, double mu, double s, double g) {
  const double b = d7::sqrt_cr(2 / M_PI);
  const double sg = (g >= 0) ? 1.0 : -1.0;
  const double cc = sg * R_pow(2 * (sg * g) / (4 - M_PI), 1.0 / 3.0);
  const double max_skew = (4 - M_PI) / 2 * R_pow(b / d7::sqrt_cr(1 - b * b), 3.0);
  const double omd2 = -std::expm1((2.0 / 3.0) * std::log(sg * g / max_skew));
  const double xi = mu - s * cc;
  const double omega = s * d7::sqrt_cr(1 + cc * cc);
  const double alpha = cc / (b * d7::sqrt_cr(omd2));
  const double z = (y - xi) / omega;
  return std::log(2.0) - std::log(omega) + R::dnorm4(z, 0.0, 1.0, 1) +
    R::pnorm5(alpha * z, 0.0, 1.0, 1, 1);
}

// the diagonal of the expected information and of its derivative in the
// same parameter for |g| < SN2_GE, from the series of
// skewnormal2_expected_series_cpp and skewnormal2_dexpected1_series_cpp
inline double skewnormal2_series_info(int k, double s, double r) {
  return sn2_spow(s, SN2_EP[k], 0) * sn2_eser(SN2_E[k], r, 0);
}

inline double skewnormal2_series_dinfo(int k, double s, double r) {
  if (k == 0) return 0.0;
  if (k == 1) return sn2_spow(s, SN2_EP[k], 1) * sn2_eser(SN2_E[k], r, 0);
  return sn2_spow(s, SN2_EP[k], 0) * sn2_eser(SN2_E[k], r, 1);
}

// the routers of the scalar registry: th = (mu, sigma, gamma1)
inline void skewnormal2_score_curv(int k, double y, const double* th,
                                   double* out) {
  const double m = th[0], s = th[1], gv = th[2];
  const Sn2ShapeAll P = sn2_shape_all(gv);
  const double r = P.r, Dq = P.Dq;
  const double w = (y - m) / s;
  const double X = r / Dq * (w + r) / d7::sqrt_cr(1.0 + r * r);
  const double Z1 = sn2_zeta1(X);
  if (k == 0) {
    out[0] = skewnormal2_score_mu(sn2_grad_G1_0(w, X, Z1, P), s);
    out[1] = skewnormal2_hess_mu_mu(sn2_G2_0(w, X, Z1, P), s);
  } else if (k == 1) {
    out[0] = skewnormal2_score_sigma(sn2_grad_G1_0(w, X, Z1, P), w, s);
    out[1] = skewnormal2_hess_sigma_sigma(sn2_G1_0(w, X, Z1, P),
                                          sn2_G2_0(w, X, Z1, P), w, s);
  } else {
    const bool ser = sn2_series_region(X, r);
    out[0] = skewnormal2_score_gamma1(
      sn2_G0_1(w, X, Z1, P, ser ? sn2_G_series(0, 1, w, r) : 0.0));
    out[1] = skewnormal2_hess_gamma1_gamma1(
      sn2_G0_2(w, X, Z1, P, ser ? sn2_G_series(0, 2, w, r) : 0.0));
  }
}

// the standardized diagonal pair of parameter k at shape = (gamma1), by the
// rule of pt_loc_scale.h
inline void skewnormal2_quad_diag(int k, const double* shape, double* out,
                                  bool want_d) {
  const double gv = shape[0];
  const Sn2ShapeAll P = sn2_shape_all(gv);
  const double r = P.r, Dq = P.Dq;
  const LocScaleRule R = loc_scale_rule();
  LocScaleSum S;
  for (int j = 0; j < R.n; ++j) {
    const double w = R.x[j];
    const double fw = std::exp(skewnormal2_logpdf(w, 0.0, 1.0, gv)) * R.w[j];
    if (fw == 0.0) continue;
    const double X = r / Dq * (w + r) / d7::sqrt_cr(1.0 + r * r);
    const double Z1 = sn2_zeta1(X);
    double g, H, T3 = 0.0;
    if (k == 0) {
      const double G2_0 = sn2_G2_0(w, X, Z1, P);
      H = skewnormal2_hess_mu_mu(G2_0, 1.0);
      if (!want_d) { S.add0(fw, H); continue; }
      g = skewnormal2_score_mu(sn2_grad_G1_0(w, X, Z1, P), 1.0);
      T3 = skewnormal2_d3_mu_mu_mu(sn2_G3_0(w, X, Z1, P), 1.0);
    } else if (k == 1) {
      const double G1_0 = sn2_G1_0(w, X, Z1, P);
      const double G2_0 = sn2_G2_0(w, X, Z1, P);
      H = skewnormal2_hess_sigma_sigma(G1_0, G2_0, w, 1.0);
      if (!want_d) { S.add0(fw, H); continue; }
      g = skewnormal2_score_sigma(sn2_grad_G1_0(w, X, Z1, P), w, 1.0);
      T3 = skewnormal2_d3_sigma_sigma_sigma(G1_0, G2_0,
                                            sn2_G3_0(w, X, Z1, P), w, 1.0);
    } else {
      const bool ser = sn2_series_region(X, r);
      H = skewnormal2_hess_gamma1_gamma1(
        sn2_G0_2(w, X, Z1, P, ser ? sn2_G_series(0, 2, w, r) : 0.0));
      if (!want_d) { S.add0(fw, H); continue; }
      g = skewnormal2_score_gamma1(
        sn2_G0_1(w, X, Z1, P, ser ? sn2_G_series(0, 1, w, r) : 0.0));
      T3 = skewnormal2_d3_gamma1_gamma1_gamma1(
        sn2_G0_3(w, X, Z1, P, ser ? sn2_G_series(0, 3, w, r) : 0.0));
    }
    S.add(fw, g, H, T3);
  }
  S.result(out);
}

inline void skewnormal2_info_dinfo(int k, double y, const double* th,
                                   double* out) {
  (void) y;
  const double s = th[1], gv = th[2];
  if (std::fabs(gv) < SN2_GE) {
    const double r = sn2_r(gv);
    out[0] = skewnormal2_series_info(k, s, r);
    out[1] = skewnormal2_series_dinfo(k, s, r);
    return;
  }
  double std2[2];
  loc_scale_cached(4, k, th + 2, 1, [&](double* v) {
    skewnormal2_quad_diag(k, th + 2, v, true);
  }, std2);
  loc_scale_unstandardize(k, s, std2, out);
}

} // namespace d7

#endif
