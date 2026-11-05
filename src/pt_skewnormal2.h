#ifndef D7_PT_SKEWNORMAL2_H
#define D7_PT_SKEWNORMAL2_H

#include <Rcpp.h>
#include <cmath>
#include "pt_loc_scale.h"
#include "pt_sqrt.h"
#include "pt_sn2_fdb.h"
#include "pt_parts.h"

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

// the same temporaries from r and Dq already formed (the parts route)
inline Sn2ShapeGrad sn2_shape_grad_rd(double r, double Dq) {
  Sn2ShapeGrad P;
  P.r = r;
  P.Dq = Dq;
  P.p0 = std::pow(P.r, 2) + 1;
  P.p1 = std::pow(P.r, 2);
  P.p2 = P.p1 + 1;
  P.p3 = 1.0/P.p2;
  P.p4 = 1.0/P.Dq;
  P.p5 = std::pow(P.p2, -1.0/2.0);
  return P;
}

inline Sn2ShapeGrad sn2_shape_grad(double g) {
  return sn2_shape_grad_rd(sn2_r(g), sn2_Dq(g));
}

struct Sn2ShapeHess {
  double r, Dq;
  double p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17;
};

inline Sn2ShapeHess sn2_shape_hess_rd(double r, double Dq) {
  Sn2ShapeHess P;
  P.r = r;
  P.Dq = Dq;
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

inline Sn2ShapeHess sn2_shape_hess(double g) {
  return sn2_shape_hess_rd(sn2_r(g), sn2_Dq(g));
}

struct Sn2ShapeD3 {
  double r, Dq;
  double p1, p2, p3, p4, p5, p6, p7, p8, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37;
};

inline Sn2ShapeD3 sn2_shape_d3_rd(double r, double Dq) {
  Sn2ShapeD3 P;
  P.r = r;
  P.Dq = Dq;
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

inline Sn2ShapeD3 sn2_shape_d3(double g) {
  return sn2_shape_d3_rd(sn2_r(g), sn2_Dq(g));
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

// every shape-only temporary of the sets above (p0 to p37, with the same
// expressions), for the full quadrature, which reads all the G_ak at once
struct Sn2ShapeFull {
  double r, Dq;
  double p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37;
};

inline Sn2ShapeFull sn2_shape_full(double g) {
  Sn2ShapeFull P;
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

template <class P_t>
inline double sn2_grad_G1_0(double w, double X, double Z1, const P_t& P) {
  // far left of the step the forms below cancel; there by composition
  if (X < -1.0 && !sn2_series_region(X, P.r)) {
    double GG[6][6];
    sn2_G_fdb(w, P.r, P.Dq, 1, GG);
    return GG[1][0];
  }
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
  // far left of the step the forms below cancel; there by composition
  if (X < -1.0 && !sn2_series_region(X, P.r)) {
    double GG[6][6];
    sn2_G_fdb(w, P.r, P.Dq, 1, GG);
    return GG[1][0];
  }
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
  // far left of the step the forms below cancel; there by composition
  if (X < -1.0 && !sn2_series_region(X, P.r)) {
    double GG[6][6];
    sn2_G_fdb(w, P.r, P.Dq, 2, GG);
    return GG[2][0];
  }
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
  // far left of the step the forms below cancel; there by composition
  if (X < -1.0 && !sn2_series_region(X, P.r)) {
    double GG[6][6];
    sn2_G_fdb(w, P.r, P.Dq, 3, GG);
    return GG[3][0];
  }
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
  // far left of the step the forms below cancel; there by composition
  if (X < -1.0 && !sn2_series_region(X, P.r)) {
    double GG[6][6];
    sn2_G_fdb(w, P.r, P.Dq, 1, GG);
    return GG[0][1];
  }
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
  // far left of the step the forms below cancel; there by composition
  if (X < -1.0 && !sn2_series_region(X, P.r)) {
    double GG[6][6];
    sn2_G_fdb(w, P.r, P.Dq, 2, GG);
    return GG[0][2];
  }
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
  // far left of the step the forms below cancel; there by composition
  if (X < -1.0 && !sn2_series_region(X, P.r)) {
    double GG[6][6];
    sn2_G_fdb(w, P.r, P.Dq, 3, GG);
    return GG[0][3];
  }
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
  // far left of the step the forms below cancel; there by composition
  if (X < -1.0 && !sn2_series_region(X, P.r)) {
    double GG[6][6];
    sn2_G_fdb(w, P.r, P.Dq, 2, GG);
    return GG[1][1];
  }
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
  // far left of the step the forms below cancel; there by composition
  if (X < -1.0 && !sn2_series_region(X, P.r)) {
    double GG[6][6];
    sn2_G_fdb(w, P.r, P.Dq, 2, GG);
    return GG[1][1];
  }
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
  // far left of the step the forms below cancel; there by composition
  if (X < -1.0 && !sn2_series_region(X, P.r)) {
    double GG[6][6];
    sn2_G_fdb(w, P.r, P.Dq, 3, GG);
    return GG[2][1];
  }
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
  // far left of the step the forms below cancel; there by composition
  if (X < -1.0 && !sn2_series_region(X, P.r)) {
    double GG[6][6];
    sn2_G_fdb(w, P.r, P.Dq, 3, GG);
    return GG[1][2];
  }
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

// the entry c of the expected Hessian, in the pair order mm, ss, gg, ms, mg,
// sg, and its derivative in parameter d (0 mu, 1 sigma, 2 gamma1), for
// |g| < SN2_GE
inline double skewnormal2_series_pair_info(int c, double s, double r) {
  return sn2_spow(s, SN2_EP[c], 0) * sn2_eser(SN2_E[c], r, 0);
}

inline double skewnormal2_series_pair_dinfo(int c, int d, double s, double r) {
  if (d == 0) return 0.0;
  if (d == 1) return sn2_spow(s, SN2_EP[c], 1) * sn2_eser(SN2_E[c], r, 0);
  return sn2_spow(s, SN2_EP[c], 0) * sn2_eser(SN2_E[c], r, 1);
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
  // the rule is centered at the step of the density, w = -r, where its
  // nodes cluster; at a skewness near the ceiling the step is too sharp for
  // nodes clustered at w = 0
  for (int j = 0; j < R.n; ++j) {
    const double w = R.x[j] - r;
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

// The off-diagonal components of the observed Hessian, from
// skewnormal2_hessian_cpp (skewnormal2.cpp): u0 = s^-2 and u1 = G_11/s.
inline double skewnormal2_hess_mu_sigma(double G1_0, double G2_0, double w,
                                        double s) {
  const double u0 = std::pow(s, -2);
  return u0 * (G1_0 + G2_0 * w);
}

inline double skewnormal2_hess_mu_gamma1(double G1_1, double s) {
  const double u1 = G1_1 / s;
  return -u1;
}

inline double skewnormal2_hess_sigma_gamma1(double G1_1, double w, double s) {
  const double u1 = G1_1 / s;
  return -u1 * w;
}

// The routers of the joint registry, in (mu, sigma, gamma1). The expected
// information is a series below sn2_ge() and a quadrature above it
// (skewnormal2_info_dinfo()), so it is not exact here and only the gradient
// and the observed Hessian are written. The observed Hessian is rejected at
// gamma1 = 0 by distrib_hessian(), and the router computes what the kernel
// computes there.
inline void skewnormal2_joint_grad(double y, const double* th, double* out) {
  const double m = th[0], s = th[1], gv = th[2];
  const Sn2ShapeGrad P = sn2_shape_grad(gv);
  const double r = P.r, Dq = P.Dq;
  const double w = (y - m) / s;
  const double X = r / Dq * (w + r) / d7::sqrt_cr(1.0 + r * r);
  const double Z1 = sn2_zeta1(X);
  double G1_0 = sn2_grad_G1_0(w, X, Z1, P);
  const double gs01 = sn2_series_region(X, r) ? sn2_G_series(0, 1, w, r) : 0.0;
  double G0_1 = sn2_G0_1(w, X, Z1, P, gs01);
  // far left of the step the forms above cancel; there by composition
  if (X < -1.0 && !sn2_series_region(X, r)) {
    double GG[6][6];
    sn2_G_fdb(w, r, Dq, 1, GG);
    G1_0 = GG[1][0];
    G0_1 = GG[0][1];
  }
  out[0] = skewnormal2_score_mu(G1_0, s);
  out[1] = skewnormal2_score_sigma(G1_0, w, s);
  out[2] = skewnormal2_score_gamma1(G0_1);
}

inline void skewnormal2_joint_hess(double y, const double* th, double* out) {
  const double m = th[0], s = th[1], gv = th[2];
  const Sn2ShapeHess P = sn2_shape_hess(gv);
  const double r = P.r, Dq = P.Dq;
  const double w = (y - m) / s;
  const double X = r / Dq * (w + r) / d7::sqrt_cr(1.0 + r * r);
  const double Z1 = sn2_zeta1(X);
  const bool ser = sn2_series_region(X, r);
  const double gs11 = ser ? sn2_G_series(1, 1, w, r) : 0.0;
  const double gs02 = ser ? sn2_G_series(0, 2, w, r) : 0.0;
  double G1_1 = sn2_hess_G1_1(w, X, Z1, P, gs11);
  double G1_0 = sn2_G1_0(w, X, Z1, P);
  double G2_0 = sn2_G2_0(w, X, Z1, P);
  double G0_2 = sn2_G0_2(w, X, Z1, P, gs02);
  // far left of the step the forms above cancel; there by composition
  if (X < -1.0 && !ser) {
    double GG[6][6];
    sn2_G_fdb(w, r, Dq, 2, GG);
    G1_0 = GG[1][0];
    G2_0 = GG[2][0];
    G1_1 = GG[1][1];
    G0_2 = GG[0][2];
  }
  out[0] = skewnormal2_hess_mu_mu(G2_0, s);
  out[4] = skewnormal2_hess_sigma_sigma(G1_0, G2_0, w, s);
  out[8] = skewnormal2_hess_gamma1_gamma1(G0_2);
  out[1] = out[3] = skewnormal2_hess_mu_sigma(G1_0, G2_0, w, s);
  out[2] = out[6] = skewnormal2_hess_mu_gamma1(G1_1, s);
  out[5] = out[7] = skewnormal2_hess_sigma_gamma1(G1_1, w, s);
}

// The third derivatives with a mixed index, from skewnormal2_deriv3_cpp
// (skewnormal2.cpp)
inline double skewnormal2_d3_mu_mu_sigma(double G2_0, double G3_0, double w,
                                         double s) {
  const double u0 = std::pow(s, -3);
  return -u0*(2*G2_0 + G3_0*w);
}

inline double skewnormal2_d3_mu_mu_gamma1(double G2_1, double s) {
  const double u1 = std::pow(s, -2);
  return G2_1*u1;
}

inline double skewnormal2_d3_mu_sigma_sigma(double G1_0, double G2_0,
                                            double G3_0, double w, double s) {
  const double u0 = std::pow(s, -3);
  const double u2 = std::pow(w, 2);
  return -u0*(2*G1_0 + 4*G2_0*w + G3_0*u2);
}

inline double skewnormal2_d3_mu_sigma_gamma1(double G1_1, double G2_1,
                                             double w, double s) {
  const double u1 = std::pow(s, -2);
  const double u3 = G2_1*w;
  return u1*(G1_1 + u3);
}

inline double skewnormal2_d3_mu_gamma1_gamma1(double G1_2, double s) {
  const double u4 = G1_2/s;
  return -u4;
}

inline double skewnormal2_d3_sigma_sigma_gamma1(double G1_1, double G2_1,
                                                double w, double s) {
  const double u1 = std::pow(s, -2);
  const double u3 = G2_1*w;
  return u1*w*(2*G1_1 + u3);
}

inline double skewnormal2_d3_sigma_gamma1_gamma1(double G1_2, double w,
                                                 double s) {
  const double u4 = G1_2/s;
  return -u4*w;
}

// the ten distinct third derivatives, in the order (mmm, mms, mmg, mss, msg,
// mgg, sss, ssg, sgg, ggg), into the P = 3 array T[a + 3 b + 9 c] with every
// permutation of the indices filled
inline void skewnormal2_fill_t3(const double* v, double* T) {
  static const int idx[10][3] = {{0, 0, 0}, {0, 0, 1}, {0, 0, 2}, {0, 1, 1},
                                 {0, 1, 2}, {0, 2, 2}, {1, 1, 1}, {1, 1, 2},
                                 {1, 2, 2}, {2, 2, 2}};
  static const int perm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2},
                                 {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
  for (int q = 0; q < 10; ++q)
    for (int p = 0; p < 6; ++p)
      T[idx[q][perm[p][0]] + 3 * idx[q][perm[p][1]] +
        9 * idx[q][perm[p][2]]] = v[q];
}

// The full standardized pair at shape = (gamma1), by the rule of
// pt_loc_scale.h centered at the step w = -r, as skewnormal2_quad_diag()
// loops: at each node the score g, the Hessian H and the third derivatives
// T3 in (mu, sigma, gamma1) at location 0 and scale 1, added with the density
// weight. I is 3 x 3 and D is 3 x 3 x 3, column-major. The score takes G_10
// in its own form (sn2_grad_G1_0), the Hessian and the third derivatives take
// the forms of their kernels; far left of the step, where those cancel, all
// the G_ak come from one composition of order three.
inline void skewnormal2_quad_full(const double* shape, double* I, double* D) {
  const double gv = shape[0];
  const Sn2ShapeFull P = sn2_shape_full(gv);
  const double r = P.r, Dq = P.Dq;
  const LocScaleRule R = loc_scale_rule();
  LocScaleFull S(3);
  for (int j = 0; j < R.n; ++j) {
    const double w = R.x[j] - r;
    const double fw = std::exp(skewnormal2_logpdf(w, 0.0, 1.0, gv)) * R.w[j];
    if (fw == 0.0) continue;
    const double X = r / Dq * (w + r) / d7::sqrt_cr(1.0 + r * r);
    const double Z1 = sn2_zeta1(X);
    const bool ser = sn2_series_region(X, r);
    double G1_0g, G1_0, G2_0, G3_0, G0_1, G0_2, G0_3, G1_1h, G1_1, G2_1, G1_2;
    if (X < -1.0 && !ser) {
      double GG[6][6];
      sn2_G_fdb(w, r, Dq, 3, GG);
      G1_0g = G1_0 = GG[1][0];
      G2_0 = GG[2][0];
      G3_0 = GG[3][0];
      G0_1 = GG[0][1];
      G0_2 = GG[0][2];
      G0_3 = GG[0][3];
      G1_1h = G1_1 = GG[1][1];
      G2_1 = GG[2][1];
      G1_2 = GG[1][2];
    } else {
      G1_0g = sn2_grad_G1_0(w, X, Z1, P);
      G1_0 = sn2_G1_0(w, X, Z1, P);
      G2_0 = sn2_G2_0(w, X, Z1, P);
      G3_0 = sn2_G3_0(w, X, Z1, P);
      G0_1 = sn2_G0_1(w, X, Z1, P, ser ? sn2_G_series(0, 1, w, r) : 0.0);
      G0_2 = sn2_G0_2(w, X, Z1, P, ser ? sn2_G_series(0, 2, w, r) : 0.0);
      G0_3 = sn2_G0_3(w, X, Z1, P, ser ? sn2_G_series(0, 3, w, r) : 0.0);
      G1_1h = sn2_hess_G1_1(w, X, Z1, P, ser ? sn2_G_series(1, 1, w, r) : 0.0);
      G1_1 = sn2_G1_1(w, X, Z1, P, ser ? sn2_G_series(1, 1, w, r) : 0.0);
      G2_1 = sn2_G2_1(w, X, Z1, P, ser ? sn2_G_series(2, 1, w, r) : 0.0);
      G1_2 = sn2_G1_2(w, X, Z1, P, ser ? sn2_G_series(1, 2, w, r) : 0.0);
    }
    double g[3], H[9], T3[27];
    g[0] = skewnormal2_score_mu(G1_0g, 1.0);
    g[1] = skewnormal2_score_sigma(G1_0g, w, 1.0);
    g[2] = skewnormal2_score_gamma1(G0_1);
    H[0] = skewnormal2_hess_mu_mu(G2_0, 1.0);
    H[4] = skewnormal2_hess_sigma_sigma(G1_0, G2_0, w, 1.0);
    H[8] = skewnormal2_hess_gamma1_gamma1(G0_2);
    H[1] = H[3] = skewnormal2_hess_mu_sigma(G1_0, G2_0, w, 1.0);
    H[2] = H[6] = skewnormal2_hess_mu_gamma1(G1_1h, 1.0);
    H[5] = H[7] = skewnormal2_hess_sigma_gamma1(G1_1h, w, 1.0);
    double v[10];
    v[0] = skewnormal2_d3_mu_mu_mu(G3_0, 1.0);
    v[1] = skewnormal2_d3_mu_mu_sigma(G2_0, G3_0, w, 1.0);
    v[2] = skewnormal2_d3_mu_mu_gamma1(G2_1, 1.0);
    v[3] = skewnormal2_d3_mu_sigma_sigma(G1_0, G2_0, G3_0, w, 1.0);
    v[4] = skewnormal2_d3_mu_sigma_gamma1(G1_1, G2_1, w, 1.0);
    v[5] = skewnormal2_d3_mu_gamma1_gamma1(G1_2, 1.0);
    v[6] = skewnormal2_d3_sigma_sigma_sigma(G1_0, G2_0, G3_0, w, 1.0);
    v[7] = skewnormal2_d3_sigma_sigma_gamma1(G1_1, G2_1, w, 1.0);
    v[8] = skewnormal2_d3_sigma_gamma1_gamma1(G1_2, w, 1.0);
    v[9] = skewnormal2_d3_gamma1_gamma1_gamma1(G0_3);
    skewnormal2_fill_t3(v, T3);
    S.add(fw, g, H, T3);
  }
  S.result(I, D);
}

// The routers of the joint registry for the expected Hessian and its
// derivative, in (mu, sigma, gamma1): P x P and P x P x P, column-major, both
// triangles filled. For |g| < SN2_GE they are the series of
// skewnormal2_expected_series_cpp and skewnormal2_dexpected1_series_cpp (the
// pairs in the order mm, ss, gg, ms, mg, sg; the information does not depend
// on mu), and above it the full quadrature, kept per thread for the last
// shape (tag 4 of loc_scale_diag_cpp()).
inline void skewnormal2_joint_info(double y, const double* th, double* out) {
  (void) y;
  const double s = th[1], gv = th[2];
  if (std::fabs(gv) < SN2_GE) {
    static const int ia[6] = {0, 1, 2, 0, 0, 1}, ib[6] = {0, 1, 2, 1, 2, 2};
    const double r = sn2_r(gv);
    for (int c = 0; c < 6; ++c)
      out[ia[c] + 3 * ib[c]] = out[ib[c] + 3 * ia[c]] =
        skewnormal2_series_pair_info(c, s, r);
    return;
  }
  const LocScaleFullSlot& slot = loc_scale_full_cached(
    4, th + 2, 1, [&](double* I, double* D) {
      skewnormal2_quad_full(th + 2, I, D);
    });
  loc_scale_full_info(slot, 3, s, out);
}

inline void skewnormal2_joint_dinfo(double y, const double* th, double* out) {
  (void) y;
  const double s = th[1], gv = th[2];
  if (std::fabs(gv) < SN2_GE) {
    static const int ia[6] = {0, 1, 2, 0, 0, 1}, ib[6] = {0, 1, 2, 1, 2, 2};
    const double r = sn2_r(gv);
    for (int c = 0; c < 6; ++c)
      for (int d = 0; d < 3; ++d)
        out[ia[c] + 3 * ib[c] + 9 * d] = out[ib[c] + 3 * ia[c] + 9 * d] =
          skewnormal2_series_pair_dinfo(c, d, s, r);
    return;
  }
  const LocScaleFullSlot& slot = loc_scale_full_cached(
    4, th + 2, 1, [&](double* I, double* D) {
      skewnormal2_quad_full(th + 2, I, D);
    });
  loc_scale_full_dinfo(slot, 3, s, out);
}

// The parts route (pt_parts.h), P = 3, (mu, sigma, gamma1). The buffer:
// [0] level, [1] y, [2] mu, [3] sigma, [4] gamma1, then, filled from level 1
// and read by every order, [5] r = sn2_r(gamma1), [6] Dq = sn2_Dq(gamma1),
// [7] X = r (w + r) / (Dq sqrt(1 + r^2)), [8] Z1 = phi(X)/Phi(X). The
// G_ak, which differ in form from one order to the next, are formed inside
// each order's function from these. The expected parts hold y as NaN and
// the parameters only.
inline void skewnormal2_parts(double y, const double* th, int level,
                              double* b) {
  const double m = th[0], s = th[1], gv = th[2];
  b[0] = level; b[1] = y; b[2] = m; b[3] = s; b[4] = gv;
  b[5] = b[6] = b[7] = b[8] = R_NaN;
  if (level >= 1) {
    const double r = sn2_r(gv);
    const double Dq = sn2_Dq(gv);
    const double w = (y - m) / s;
    const double X = r / Dq * (w + r) / d7::sqrt_cr(1.0 + r * r);
    b[5] = r; b[6] = Dq; b[7] = X; b[8] = sn2_zeta1(X);
  }
}

inline void skewnormal2_parts_expected(const double* th, int level,
                                       double* b) {
  (void) level;
  b[0] = level; b[1] = R_NaN; b[2] = th[0]; b[3] = th[1]; b[4] = th[2];
  b[5] = b[6] = b[7] = b[8] = R_NaN;
}

inline void skewnormal2_parts_grad(const double* b, double* out) {
  const double s = b[3];
  const double r = b[5], Dq = b[6], X = b[7], Z1 = b[8];
  const Sn2ShapeGrad P = sn2_shape_grad_rd(r, Dq);
  const double w = (b[1] - b[2]) / s;
  double G1_0 = sn2_grad_G1_0(w, X, Z1, P);
  const double gs01 = sn2_series_region(X, r) ? sn2_G_series(0, 1, w, r) : 0.0;
  double G0_1 = sn2_G0_1(w, X, Z1, P, gs01);
  // far left of the step the forms above cancel; there by composition
  if (X < -1.0 && !sn2_series_region(X, r)) {
    double GG[6][6];
    sn2_G_fdb(w, r, Dq, 1, GG);
    G1_0 = GG[1][0];
    G0_1 = GG[0][1];
  }
  out[0] = skewnormal2_score_mu(G1_0, s);
  out[1] = skewnormal2_score_sigma(G1_0, w, s);
  out[2] = skewnormal2_score_gamma1(G0_1);
}

inline void skewnormal2_parts_hess(const double* b, double* out) {
  const double s = b[3];
  const double r = b[5], Dq = b[6], X = b[7], Z1 = b[8];
  const Sn2ShapeHess P = sn2_shape_hess_rd(r, Dq);
  const double w = (b[1] - b[2]) / s;
  const bool ser = sn2_series_region(X, r);
  const double gs11 = ser ? sn2_G_series(1, 1, w, r) : 0.0;
  const double gs02 = ser ? sn2_G_series(0, 2, w, r) : 0.0;
  double G1_1 = sn2_hess_G1_1(w, X, Z1, P, gs11);
  double G1_0 = sn2_G1_0(w, X, Z1, P);
  double G2_0 = sn2_G2_0(w, X, Z1, P);
  double G0_2 = sn2_G0_2(w, X, Z1, P, gs02);
  // far left of the step the forms above cancel; there by composition
  if (X < -1.0 && !ser) {
    double GG[6][6];
    sn2_G_fdb(w, r, Dq, 2, GG);
    G1_0 = GG[1][0];
    G2_0 = GG[2][0];
    G1_1 = GG[1][1];
    G0_2 = GG[0][2];
  }
  out[0] = skewnormal2_hess_mu_mu(G2_0, s);
  out[4] = skewnormal2_hess_sigma_sigma(G1_0, G2_0, w, s);
  out[8] = skewnormal2_hess_gamma1_gamma1(G0_2);
  out[1] = out[3] = skewnormal2_hess_mu_sigma(G1_0, G2_0, w, s);
  out[2] = out[6] = skewnormal2_hess_mu_gamma1(G1_1, s);
  out[5] = out[7] = skewnormal2_hess_sigma_gamma1(G1_1, w, s);
}

// the expected Hessian, from the expected parts: the series below SN2_GE and
// the full quadrature above it, as skewnormal2_joint_info()
inline void skewnormal2_parts_info(const double* b, double* out) {
  const double s = b[3], gv = b[4];
  if (std::fabs(gv) < SN2_GE) {
    static const int ia[6] = {0, 1, 2, 0, 0, 1}, ib[6] = {0, 1, 2, 1, 2, 2};
    const double r = sn2_r(gv);
    for (int c = 0; c < 6; ++c)
      out[ia[c] + 3 * ib[c]] = out[ib[c] + 3 * ia[c]] =
        skewnormal2_series_pair_info(c, s, r);
    return;
  }
  const double shape[1] = {gv};
  const LocScaleFullSlot& slot = loc_scale_full_cached(
    4, shape, 1, [&](double* I, double* D) {
      skewnormal2_quad_full(shape, I, D);
    });
  loc_scale_full_info(slot, 3, s, out);
}

// the ten unique components, as skewnormal2_deriv3_cpp
inline void skewnormal2_parts_d3(const double* b, double* out) {
  const double s = b[3];
  const double r = b[5], Dq = b[6], X = b[7], Z1 = b[8];
  const Sn2ShapeD3 P = sn2_shape_d3_rd(r, Dq);
  const double w = (b[1] - b[2]) / s;
  const bool ser = sn2_series_region(X, r);
  double G1_0 = 0.0;
  double G2_0 = 0.0;
  double G3_0 = 0.0;
  double G1_1 = 0.0;
  double G2_1 = 0.0;
  double G1_2 = 0.0;
  double G0_3 = 0.0;
  G1_1 = sn2_G1_1(w, X, Z1, P, ser ? sn2_G_series(1, 1, w, r) : 0.0);
  G2_1 = sn2_G2_1(w, X, Z1, P, ser ? sn2_G_series(2, 1, w, r) : 0.0);
  G1_2 = sn2_G1_2(w, X, Z1, P, ser ? sn2_G_series(1, 2, w, r) : 0.0);
  G1_0 = sn2_G1_0(w, X, Z1, P);
  G2_0 = sn2_G2_0(w, X, Z1, P);
  G3_0 = sn2_G3_0(w, X, Z1, P);
  G0_3 = sn2_G0_3(w, X, Z1, P, ser ? sn2_G_series(0, 3, w, r) : 0.0);
  // far left of the step the forms above cancel; there by composition
  if (X < -1.0 && !ser) {
    double GG[6][6];
    sn2_G_fdb(w, r, Dq, 3, GG);
    G1_0 = GG[1][0];
    G2_0 = GG[2][0];
    G3_0 = GG[3][0];
    G1_1 = GG[1][1];
    G2_1 = GG[2][1];
    G1_2 = GG[1][2];
    G0_3 = GG[0][3];
  }
  double u[10];
  u[0] = skewnormal2_d3_mu_mu_mu(G3_0, s);
  u[1] = skewnormal2_d3_mu_mu_sigma(G2_0, G3_0, w, s);
  u[2] = skewnormal2_d3_mu_mu_gamma1(G2_1, s);
  u[3] = skewnormal2_d3_mu_sigma_sigma(G1_0, G2_0, G3_0, w, s);
  u[4] = skewnormal2_d3_mu_sigma_gamma1(G1_1, G2_1, w, s);
  u[5] = skewnormal2_d3_mu_gamma1_gamma1(G1_2, s);
  u[6] = skewnormal2_d3_sigma_sigma_sigma(G1_0, G2_0, G3_0, w, s);
  u[7] = skewnormal2_d3_sigma_sigma_gamma1(G1_1, G2_1, w, s);
  u[8] = skewnormal2_d3_sigma_gamma1_gamma1(G1_2, w, s);
  u[9] = skewnormal2_d3_gamma1_gamma1_gamma1(G0_3);
  sym_fill(3, 3, u, out);
}

// the fifteen unique components, as skewnormal2_deriv4_cpp: the G_ak of
// order four are its generated expressions, with sn2_series() read through
// sn2_G_series() (the same sums to the last bit)
inline void skewnormal2_parts_d4(const double* b, double* out) {
  const double s = b[3];
  const double r = b[5], Dq = b[6], X = b[7], Z1 = b[8];
  const double w = (b[1] - b[2]) / s;
  double G1_0 = 0.0;
  double G2_0 = 0.0;
  double G3_0 = 0.0;
  double G4_0 = 0.0;
  double G1_1 = 0.0;
  double G2_1 = 0.0;
  double G3_1 = 0.0;
  double G1_2 = 0.0;
  double G2_2 = 0.0;
  double G1_3 = 0.0;
  double G0_4 = 0.0;
  if (std::fabs(X) < SN2_XC && std::fabs(r) < SN2_RC) {
    G1_1 = sn2_G_series(1, 1, w, r);
    G2_1 = sn2_G_series(2, 1, w, r);
    G3_1 = sn2_G_series(3, 1, w, r);
    G1_2 = sn2_G_series(1, 2, w, r);
    G2_2 = sn2_G_series(2, 2, w, r);
    G1_3 = sn2_G_series(1, 3, w, r);
    G0_4 = sn2_G_series(0, 4, w, r);
    const double e0 = std::pow(r, 2);
    const double e1 = e0 + 1;
    const double e2 = 1.0/e1;
    const double e3 = std::pow(X, 2);
    const double e4 = 3*X;
    const double e5 = std::pow(Z1, 2);
    G1_0 = -e2*(r + w) + Z1*r/(Dq*d7::sqrt_cr(e1));
    G2_0 = -e2*(1 + Z1*e0*(X + Z1)/std::pow(Dq, 2));
    G3_0 = Z1*std::pow(r, 3)*(Z1*e4 + e3 + 2*e5 - 1)/(std::pow(Dq, 3)*std::pow(e1, 3.0/2.0));
    G4_0 = -Z1*std::pow(r, 4)*(std::pow(X, 3) + 12*X*e5 + 6*std::pow(Z1, 3) + 7*Z1*e3 - 4*Z1 - e4)/(std::pow(Dq, 4)*std::pow(e1, 2));
  } else {
    const double f0 = r + w;
    const double f1 = std::pow(r, 2);
    const double f2 = f1 + 1;
    const double f3 = 1.0/f2;
    const double f4 = f0*f3;
    const double f5 = 1.0/Dq;
    const double f6 = std::pow(f2, -1.0/2.0);
    const double f7 = std::pow(Dq, -2);
    const double f8 = f1*f7;
    const double f9 = X + Z1;
    const double f10 = Z1*f9;
    const double f11 = std::pow(r, 3);
    const double f12 = std::pow(f2, -3.0/2.0);
    const double f13 = std::pow(X, 2);
    const double f14 = 3*X;
    const double f15 = std::pow(Z1, 2);
    const double f16 = Z1*f14 + f13 + 2*f15 - 1;
    const double f17 = f12*f16;
    const double f18 = std::pow(Dq, -3);
    const double f19 = Z1*f18;
    const double f20 = f17*f19;
    const double f21 = std::pow(f2, -2);
    const double f22 = std::pow(r, 4);
    const double f23 = f21*f22;
    const double f24 = 4*Z1;
    const double f25 = std::pow(X, 3) + 12*X*f15 + 6*std::pow(Z1, 3) + 7*Z1*f13 - f14 - f24;
    const double f26 = Z1*f25/std::pow(Dq, 4);
    const double f27 = 1.0/r;
    const double f28 = 1 - 2/M_PI;
    const double f29 = f28*f8;
    const double f30 = f0*f29;
    const double f31 = 2*r;
    const double f32 = -f0*f1*f3 + f31 + w;
    const double f33 = -f30 - f32;
    const double f34 = f0*f31;
    const double f35 = f1*f3;
    const double f36 = -f28;
    const double f37 = f36*f8;
    const double f38 = f35 + f37 - 1;
    const double f39 = Z1*f6;
    const double f40 = f39*f5;
    const double f41 = -f21*f34 + f3 + f38*f40;
    const double f42 = 4 - M_PI;
    const double f43 = 1.0/f42;
    const double f44 = (2.0/3.0)*f43;
    const double f45 = 2*f21;
    const double f46 = -f1*f3 + 1;
    const double f47 = -f29 - f46;
    const double f48 = f27*f47;
    const double f49 = 2*f3;
    const double f50 = f10*f7;
    const double f51 = f49*f50;
    const double f52 = f0*f37;
    const double f53 = f1*f4 - f31 - w;
    const double f54 = f52 + f53;
    const double f55 = 1.0/f1;
    const double f56 = std::pow(f33, 2);
    const double f57 = f27*f56;
    const double f58 = 2*f55;
    const double f59 = 3*f30 + f32;
    const double f60 = -f59;
    const double f61 = f28*f7;
    const double f62 = f60*f61;
    const double f63 = f36*f7;
    const double f64 = f11*f4*f63;
    const double f65 = f11*f21;
    const double f66 = 3*f65;
    const double f67 = f0*f66;
    const double f68 = 3*f4;
    const double f69 = 2*f35 - f67 + f68*r - 2;
    const double f70 = f34*f63 + f37 - f64 + f69;
    const double f71 = f27*f70;
    const double f72 = -f54*f58 + f62 + f71;
    const double f73 = f3*f50;
    const double f74 = f72*f73;
    const double f75 = 3*f29;
    const double f76 = f46 + f75;
    const double f77 = f18*f28*f39;
    const double f78 = f27*f54;
    const double f79 = f21*r;
    const double f80 = std::pow(f2, -3);
    const double f81 = f1*f80;
    const double f82 = 3*f3;
    const double f83 = f1*f21;
    const double f84 = 2*f61;
    const double f85 = f35*f61;
    const double f86 = f84 - f85;
    const double f87 = -f82 + 3*f83 + f86;
    const double f88 = f0*f45 - 8*f0*f81 + f40*f87*r + 4*f79;
    const double f89 = -f27*f88 + f38*f51*f78 - f41*f58 - f76*f77;
    const double f90 = (4.0/9.0)/std::pow(f42, 2);
    const double f91 = 1.0/f11;
    const double f92 = f21*f26;
    const double f93 = f61*f76;
    const double f94 = f47*f58 + f87 + f93;
    const double f95 = f27*f51;
    const double f96 = std::pow(f42, -3);
    const double f97 = f27*f33;
    const double f98 = -3*f0*f3*r - 2*f1*f3 + f67 + 2;
    const double f99 = -f0*f11*f28*f3*f7 + f29 + f34*f61 + f98;
    const double f100 = 2*f27;
    const double f101 = -f72;
    const double f102 = f58*f60;
    const double f103 = 6*f52 + f53;
    const double f104 = f103*f84;
    const double f105 = f0*r;
    const double f106 = 3*f0*f11*f28*f3*f7 - 6*f105*f61 - f75 - f98;
    const double f107 = f106*f27;
    const double f108 = -f102 + f104 + f107 + f62;
    const double f109 = 6*f91;
    const double f110 = f106*f61;
    const double f111 = f3*r;
    const double f112 = 9*f65;
    const double f113 = 4*r;
    const double f114 = 2*f63;
    const double f115 = 18*f83;
    const double f116 = f22*f80;
    const double f117 = 15*f116;
    const double f118 = f11*f49;
    const double f119 = 5*f4;
    const double f120 = 3*f23;
    const double f121 = f120*f63;
    const double f122 = f0*f114 - f0*f115 + f0*f117 + f0*f121 + 9*f111 - f112 + f113*f63 - f118*f63 - f119*f37 + f68;
    const double f123 = f122*f27;
    const double f124 = f27*(-f109*f54 - f110 - f123 + f27*f28*f60*f7 + 4*f55*f70);
    const double f125 = f101*f58 + f108*f61 - f124;
    const double f126 = 1.0/f22;
    const double f127 = -f3 + f83 + f86;
    const double f128 = f20*f55;
    const double f129 = 6*f21;
    const double f130 = 24*f81;
    const double f131 = f105*f80;
    const double f132 = std::pow(f2, -4);
    const double f133 = f0*f11*f132;
    const double f134 = 3*f6;
    const double f135 = f5*f9;
    const double f136 = f135*f49;
    const double f137 = f19*f28;
    const double f138 = std::pow(r, 5);
    const double f139 = 4*f55;
    const double f140 = std::pow(f54, 2);
    const double f141 = f128*f140;
    const double f142 = std::pow(f54, 3);
    const double f143 = f136*f60;
    const double f144 = f104*f6;
    const double f145 = f144 + f6*f62;
    const double f146 = -f102*f6 + f107*f6 + f143*f78 + f145;
    const double f147 = f137*f146;
    const double f148 = f60*f77;
    const double f149 = f140*f73;
    const double f150 = std::pow(f0, 2);
    const double f151 = -f111 + f150*f21*r - f4 - f40*f54;
    const double f152 = 4*f0*f21*r + 2*f1*f21 + f150*f21 - 4*f150*f81 - f40*f70 - f49;
    const double f153 = -f148 - f149*f27 - f151*f58 + f152*f27;
    const double f154 = f140*f55;
    const double f155 = f51*f54;
    const double f156 = f11*f80;
    const double f157 = f150*f80;
    const double f158 = -f0*f130 + 6*f0*f21 + 24*f11*f132*f150 - f122*f40 - 8*f156 - 12*f157*r + 12*f21*r;
    const double f159 = -f106*f77 + f109*f151 - f139*f152 + f148*f27 + f154*f51 - f155*f71 + f158*f27;
    const double f160 = f135*f3;
    const double f161 = f160*f27;
    const double f162 = f6*(12*f105*f63 + 6*f37 - 6*f64 + f69);
    const double f163 = 4*f3;
    const double f164 = f54*f55;
    const double f165 = f113*f61;
    const double f166 = f120*f61;
    const double f167 = f134*(-f0*f1*f129 + 5*f0*f116 - f0*f166 - f0*f84 + f118*f61 + f119*f29 - f165 + f4 - f66 + f82*r);
    const double f168 = -f106*f139*f6 + f107*f136*f54 + f109*f6*f60 + f110*f6 - f135*f163*f164*f60 + f143*f71 - f144*f27 + f162*f84 + f167*f27 - f27*f6*f62;
    const double f169 = f0*f138;
    G1_0 = Z1*f5*f6*r - f4;
    G2_0 = -f3*(f10*f8 + 1);
    G3_0 = f11*f20;
    G4_0 = -f23*f26;
    G1_1 = f27*f44*(Z1*f3*f33*f7*f9 - f27*f41);
    G2_1 = f44*(-f20*f33 + f27*f45 + f48*f51);
    G3_1 = 2*f19*f43*(-f17*f38 + (1.0/3.0)*f21*f25*f5*f54*r);
    G1_2 = f55*f90*(f20*f57 - f27*f89 + f74);
    G2_2 = f27*f90*(2*Z1*f12*f16*f18*f33*f47*f55 + Z1*f12*f16*f18*(2*f33*f47*f55 - f72) - f49*f91*(f3*(4*f35 - 1) + std::pow(f38, 2)*f50 + f49) - f57*f92 - f94*f95);
    G1_3 = (8.0/27.0)*f91*f96*(f125*f73 + f20*f97*(f100*(-f97 - f99) - f33*f58 - f59*f84 + f72) + f27*(3*Z1*f127*f18*f28*f6 + 2*Z1*f27*f3*f33*f7*f87*f9 + Z1*f27*f3*f33*f7*f9*f94 + 4*Z1*f3*f33*f47*f7*f9*f91 - 6*f126*f41 - 3*f128*f47*f56 - f137*(-f127*f134 - f136*f48*f59 - f136*f76*f97 + 2*f28*f6*f7*(-6*f29 - f46) + 2*f55*f6*f76 - f6*f93) + f47*f51*f55*f99 - f48*f74 - f55*f76*f77 + 2*f55*f89 - f55*(-f129 + f130 + 24*f131 - 48*f133 + f40*(f114 - f115 + f117 + f121 - 5*f35*f63 + f82)) - 4*f88*f91) + std::pow(f33, 3)*f55*f92);
    G0_4 = -8.0/81.0*f96*(f125*f73*f78 + f137*(-f101*f161*f60 + f108*f160*f78 - f146*f58 + f146*f61 + 3*f154*f17*f60*f7 + f168*f27 - f84*(-f103*f136*f78 + 2*f103*f55*f6 - f134*f61*(-10*f30 - f32) - f145 - f161*std::pow(f60, 2) - f162*f27)) + f141*(f100*(2*f0*f36*f7*r + f1*f36*f7 - f64 - f78 - f98) - f139*f54 + 3*f62 + f71) - f27*(-f101*f164*f51 + f101*f71*f73 + f109*f153 + f124*f54*f73 - f137*f168 - f139*f159 - 3*f141*f70 + f142*f17*f18*f24*f91 + f147*f27 + f27*(2*Z1*f106*f18*f27*f28*f6 + 8*Z1*f3*f54*f55*f7*f70*f9 - f102*f77 - f109*f149 - f123*f155 - 24*f126*f151 - f137*f167 + 18*f152*f91 - 6*f158*f55 + 3*f27*(48*f1*f132*f150 + f129 - 32*f131 + 16*f132*f22 + 64*f133 - 64*f150*f22/std::pow(f2, 5) - 4*f157 + f40*(f0*f112*f61 - 50*f0*f156 + 15*f0*f79 - 20*f116 + 35*f132*f169 - f163 - f165*f4 + f166 - 5*f169*f61*f80 + 24*f83 + f84 - 5*f85) - 32*f81) - std::pow(f70, 2)*f95)) + std::pow(f54, 4)*f91*f92 + f58*(Z1*f101*f27*f3*f54*f7*f9 - f128*f142 - f147 - f153*f58 + f159*f27))/(f138*(2 - 1.0/2.0*M_PI));
  }
  // far left of the step the forms above cancel; there by composition
  if (X < -1.0 && !sn2_series_region(X, r)) {
    double GG[6][6];
    sn2_G_fdb(w, r, Dq, 4, GG);
    G1_0 = GG[1][0];
    G2_0 = GG[2][0];
    G3_0 = GG[3][0];
    G4_0 = GG[4][0];
    G1_1 = GG[1][1];
    G2_1 = GG[2][1];
    G3_1 = GG[3][1];
    G1_2 = GG[1][2];
    G2_2 = GG[2][2];
    G1_3 = GG[1][3];
    G0_4 = GG[0][4];
  }
  const double u0 = std::pow(s, -4);
  const double u1 = std::pow(s, -3);
  const double u2 = 6*w;
  const double u3 = std::pow(w, 2);
  const double u4 = std::pow(s, -2);
  const double u5 = std::pow(w, 3);
  const double u6 = G3_1*u3;
  const double u7 = G2_2*w;
  const double u8 = G1_3/s;
  double v[15];
  v[0] = G4_0*u0;
  v[1] = u0*(3*G3_0 + G4_0*w);
  v[2] = -G3_1*u1;
  v[3] = u0*(6*G2_0 + G3_0*u2 + G4_0*u3);
  v[4] = -u1*(2*G2_1 + G3_1*w);
  v[5] = G2_2*u4;
  v[6] = u0*(6*G1_0 + 18*G2_0*w + 9*G3_0*u3 + G4_0*u5);
  v[7] = -u1*(2*G1_1 + 4*G2_1*w + u6);
  v[8] = u4*(G1_2 + u7);
  v[9] = -u8;
  v[10] = u0*(24*G1_0*w + 36*G2_0*u3 + 12*G3_0*u5 + G4_0*std::pow(w, 4) + 6);
  v[11] = -u1*w*(6*G1_1 + G2_1*u2 + u6);
  v[12] = u4*w*(2*G1_2 + u7);
  v[13] = -u8*w;
  v[14] = G0_4;
  sym_fill(3, 4, v, out);
}

} // namespace d7

#endif
