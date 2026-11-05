#ifndef D7_PT_GENGAMMA2_H
#define D7_PT_GENGAMMA2_H

#include <Rcpp.h>
#include <cmath>
#include "pt_parts.h"

// The generalized gamma in its mean m, shape d and power p: one function per
// component and order for the quantities the scalar registry reads, lifted
// from the generated kernels in gengamma2.cpp (which records the
// derivation) and called by them and by the registry (d7_ccallable.cpp).
// Each component is a polynomial in the data quantities Z, R or Q, X whose
// coefficients depend on (d, p) alone; gengamma2_coef_<component>() forms
// that component's coefficients, in double or in double-double, from the
// polygamma values the caller passes, and the coefficient generators of
// the kernels call these pieces. gengamma2_coefs_<component>() computes
// the polygammas that one component needs and then its coefficients.

// --- double-double arithmetic (Dekker, Knuth; as in the QD library) --------
struct dd {
  double hi, lo;
  dd(double h = 0.0, double l = 0.0) : hi(h), lo(l) {}
};
static inline dd dd_two_sum(double a, double b) {
  const double s = a + b, bb = s - a;
  return dd(s, (a - (s - bb)) + (b - bb));
}
static inline dd dd_quick(double a, double b) {
  const double s = a + b;
  return dd(s, b - (s - a));
}
static inline dd operator+(const dd& a, const dd& b) {
  dd s = dd_two_sum(a.hi, b.hi);
  const dd t = dd_two_sum(a.lo, b.lo);
  s.lo += t.hi;
  s = dd_quick(s.hi, s.lo);
  s.lo += t.lo;
  return dd_quick(s.hi, s.lo);
}
static inline dd operator-(const dd& a) { return dd(-a.hi, -a.lo); }
static inline dd operator-(const dd& a, const dd& b) { return a + (-b); }
static inline dd operator*(const dd& a, const dd& b) {
  const double ph = a.hi * b.hi;
  double pl = std::fma(a.hi, b.hi, -ph);
  pl += a.hi * b.lo + a.lo * b.hi;
  return dd_quick(ph, pl);
}
static inline dd operator/(const dd& a, const dd& b) {
  const double q1 = a.hi / b.hi;
  dd r = a - b * dd(q1);
  const double q2 = r.hi / b.hi;
  r = r - b * dd(q2);
  const double q3 = r.hi / b.hi;
  return dd_quick(q1, q2) + dd(q3);
}
static inline dd operator+(const dd& a, double b) { return a + dd(b); }
static inline dd operator+(double a, const dd& b) { return dd(a) + b; }
static inline dd operator-(const dd& a, double b) { return a - dd(b); }
static inline dd operator-(double a, const dd& b) { return dd(a) - b; }
static inline dd operator*(const dd& a, double b) { return a * dd(b); }
static inline dd operator*(double a, const dd& b) { return dd(a) * b; }
static inline dd operator/(const dd& a, double b) { return a / dd(b); }
static inline dd operator/(double a, const dd& b) { return dd(a) / b; }
static inline dd dd_rat(double a, double b) { return dd(a) / dd(b); }
static inline dd dd_powi(dd x, int n) {
  if (n < 0) return dd(1.0) / dd_powi(x, -n);
  dd r(1.0);
  while (n) {
    if (n & 1) r = r * x;
    x = x * x;
    n >>= 1;
  }
  return r;
}

// log of a double-double: r = 2^e f, f in [0.75, 1.5), and
// log f = 2 atanh(w), w = (f - 1)/(f + 1), |w| <= 0.2
static inline dd dd_log(const dd& r) {
  int e = 0;
  const double fr = std::frexp(r.hi, &e);   // r.hi = fr 2^e, fr in [0.5, 1)
  double sc = std::ldexp(1.0, -e);
  dd f(r.hi * sc, r.lo * sc);
  if (fr < 0.75) { f = f * 2.0; e -= 1; }
  const dd w = (f - 1.0) / (f + 1.0), w2 = w * w;
  dd term = w, sum = w;
  for (int j = 1; j < 40; j++) {
    term = term * w2;
    const dd t = term / (2.0 * j + 1.0);
    sum = sum + t;
    if (std::fabs(t.hi) < 1e-34 * std::fabs(sum.hi)) break;
  }
  const dd ln2(0.6931471805599453, 2.3190468138462996e-17);
  return 2.0 * sum + ln2 * (double) e;
}
// B_2j (2j+n-1)!/(2j)!, j = 1..17, n = 1..11, as double-double
static const double GG2_BT[11][17][2] = {
  {{0.16666666666666666, 9.25185853854297e-18}, {-0.03333333333333333, -4.625929269271486e-19}, {0.023809523809523808, 1.32169407693471e-18}, {-0.03333333333333333, -4.625929269271486e-19}, {0.07575757575757576, -2.10269512239613e-18}, {-0.2531135531135531, -1.1061562736192037e-17}, {1.1666666666666667, -7.401486830834377e-17}, {-7.092156862745098, -3.274069468698501e-16}, {54.971177944862156, -1.9588897477095493e-16}, {-529.1242424242424, 6.890111377067638e-16}, {6192.123188405797, 9.226757844073186e-14}, {-86580.25311355312, 3.5926706461242705e-12}, {1425517.1666666667, -7.761021455128987e-11}, {-27298231.067816094, 1.610010519795034e-09}, {601580873.9006424, -2.6635227381825164e-08}, {-15116315767.092157, 5.011465035232843e-07}, {429614643061.1667, -2.0345052083333332e-05}},
  {{0.5, 0.0}, {-0.16666666666666666, -9.25185853854297e-18}, {0.16666666666666666, 9.25185853854297e-18}, {-0.3, -1.1102230246251566e-17}, {0.8333333333333334, -3.700743415417188e-17}, {-3.2904761904761903, -1.4380031557049647e-16}, {17.5, 0.0}, {-120.56666666666666, -3.789561257387201e-15}, {1044.452380952381, -1.0827317878249145e-14}, {-11111.60909090909, -9.921760382977398e-14}, {142418.83333333334, -9.701276818911234e-12}, {-2164506.327838828, 1.70571900112725e-11}, {38488963.5, 0.0}, {-791648700.9666667, 3.178914388020833e-08}, {18649007090.919914, -8.256920488365801e-07}, {-498838420314.0412, 2.225988051470588e-05}, {15036512507140.834, -0.0006510416666666666}},
  {{2.0, 0.0}, {-1.0, 0.0}, {1.3333333333333333, 7.401486830834377e-17}, {-3.0, 0.0}, {10.0, 0.0}, {-46.06666666666667, 3.315866100213801e-15}, {280.0, 0.0}, {-2170.2, -1.8189894035458566e-13}, {20889.04761904762, 6.929483442079453e-13}, {-244455.4, -5.820766091346741e-12}, {3418052.0, 0.0}, {-56277164.52380952, -1.419158208937872e-09}, {1077690978.0, 0.0}, {-23749461029.0, 0.0}, {596768226909.4373, -2.6422145562770563e-05}, {-16960506290677.4, 0.000390625}, {541314450257070.0, 0.0}},
  {{10.0, 0.0}, {-7.0, 0.0}, {12.0, 0.0}, {-33.0, 0.0}, {130.0, 0.0}, {-691.0, 0.0}, {4760.0, 0.0}, {-41233.8, 2.9103830456733705e-12}, {438670.0, 0.0}, {-5622474.2, 1.8626451492309571e-10}, {85451300.0, 0.0}, {-1519483442.142857, -6.811959402901785e-08}, {31253038362.0, 0.0}, {-736233291899.0, 0.0}, {19693351488011.43, -0.0011160714285714285}, {-593617720173709.0, 0.0}, {2.002863465951159e+16, -2.0}},
  {{60.0, 0.0}, {-56.0, 0.0}, {120.0, 0.0}, {-396.0, 0.0}, {1820.0, 0.0}, {-11056.0, 0.0}, {85680.0, 0.0}, {-824676.0, 0.0}, {9650740.0, 0.0}, {-134939380.8, 1.1920928955078126e-08}, {2221733800.0, 0.0}, {-42545536380.0, 0.0}, {937591150860.0, 0.0}, {-23559465340768.0, 0.0}, {669573950592388.6, -0.05357142857142857}, {-2.1370237926253524e+16, 0.0}, {7.610881170614404e+17, 36.0}},
  {{420.0, 0.0}, {-504.0, 0.0}, {1320.0, 0.0}, {-5148.0, 0.0}, {27300.0, 0.0}, {-187952.0, 0.0}, {1627920.0, 0.0}, {-17318196.0, 0.0}, {221967020.0, 0.0}, {-3373484520.0, 0.0}, {59986812600.0, 0.0}, {-1233820555020.0, 0.0}, {29065325676660.0, 0.0}, {-777462356245344.0, 0.0}, {2.34350882707336e+16, 0.0}, {-7.906988032713804e+17, -36.0}, {2.9682436565396177e+19, -516.0}},
  {{3360.0, 0.0}, {-5040.0, 0.0}, {15840.0, 0.0}, {-72072.0, 0.0}, {436800.0, 0.0}, {-3383136.0, 0.0}, {32558400.0, 0.0}, {-381000312.0, 0.0}, {5327208480.0, 0.0}, {-87710597520.0, 0.0}, {1679630752800.0, 0.0}, {-37014616650600.0, 0.0}, {930090421653120.0, 0.0}, {-2.6433720112341696e+16, 0.0}, {8.436631777464096e+17, 0.0}, {-3.0046554524312453e+19, -1624.0}, {1.187297462615847e+21, 110432.0}},
  {{30240.0, 0.0}, {-55440.0, 0.0}, {205920.0, 0.0}, {-1081080.0, 0.0}, {7425600.0, 0.0}, {-64279584.0, 0.0}, {683726400.0, 0.0}, {-8763007176.0, 0.0}, {133180212000.0, 0.0}, {-2368186133040.0, 0.0}, {48709291831200.0, 0.0}, {-1147453116168600.0, 0.0}, {3.069298391455296e+16, 0.0}, {-9.251802039319593e+17, -64.0}, {3.1215537576617157e+19, -1408.0}, {-1.1718156264481857e+21, -14184.0}, {4.867919596724973e+22, -715168.0}},
  {{302400.0, 0.0}, {-665280.0, 0.0}, {2882880.0, 0.0}, {-17297280.0, 0.0}, {133660800.0, 0.0}, {-1285591680.0, 0.0}, {15041980800.0, 0.0}, {-210312172224.0, 0.0}, {3462685512000.0, 0.0}, {-66309211725120.0, 0.0}, {1461278754936000.0, 0.0}, {-3.67184997173952e+16, 0.0}, {1.0435614530948006e+18, 0.0}, {-3.330648734155054e+19, 1792.0}, {1.1861904279114518e+21, 61184.0}, {-4.687262505792743e+22, 2578368.0}, {2.0445262306244885e+24, 104180672.0}},
  {{3326400.0, 0.0}, {-8648640.0, 0.0}, {43243200.0, 0.0}, {-294053760.0, 0.0}, {2539555200.0, 0.0}, {-26997425280.0, 0.0}, {345965558400.0, 0.0}, {-5257804305600.0, 0.0}, {93492508824000.0, 0.0}, {-1922967140028480.0, 0.0}, {4.5299641403016e+16, 0.0}, {-1.2117104906740416e+18, 0.0}, {3.652465085831802e+19, 1792.0}, {-1.23234003163737e+21, 45824.0}, {4.6261426688546625e+22, -759552.0}, {-1.9217776273750245e+24, -129167936.0}, {8.791462791685302e+25, -7599826624.0}},
  {{39916800.0, 0.0}, {-121080960.0, 0.0}, {691891200.0, 0.0}, {-5292967680.0, 0.0}, {50791104000.0, 0.0}, {-593943356160.0, 0.0}, {8303173401600.0, 0.0}, {-136702911945600.0, 0.0}, {2617790247072000.0, 0.0}, {-5.76890142008544e+16, 0.0}, {1.449588524896512e+18, 0.0}, {-4.119815668291741e+19, -1536.0}, {1.3148874308994487e+21, 130048.0}, {-4.682892120222006e+22, 2265600.0}, {1.8504570675418649e+24, 103835648.0}, {-8.071466034975104e+25, 4238623104.0}, {3.8682436283415325e+27, -59514464512.0}},
};
// B_2j/(2j), j = 1..17
static const double GG2_B0[17][2] = {{0.08333333333333333, 4.625929269271485e-18}, {-0.008333333333333333, -1.1564823173178714e-19}, {0.003968253968253968, 2.20282346155785e-19}, {-0.004166666666666667, -5.782411586589357e-20}, {0.007575757575757576, -2.1026951223961299e-19}, {-0.021092796092796094, 1.3911677399530732e-18}, {0.08333333333333333, 4.625929269271485e-18}, {-0.4432598039215686, -2.0462934179365632e-17}, {3.0539543302701198, -1.0882720820608607e-17}, {-26.456212121212122, 7.449932926454383e-16}, {281.46014492753625, -1.647635329298783e-14}, {-3607.5105463980462, -1.5347029033579816e-13}, {54827.583333333336, -2.4253192047278085e-12}, {-974936.8238505747, -4.2284185858978625e-11}, {20052695.79668808, -8.878409127275055e-10}, {-472384867.7216299, 1.5660828235102635e-08}, {12635724795.916666, 6.357828776041666e-07}};
static const double GG2_FACT[12] = {1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0, 3628800.0, 39916800.0};

// psi^(n)(z), n >= 1, z >= 40: (-1)^(n+1) [ (n-1)!/z^n + n!/(2 z^(n+1))
//   + sum_j B_2j (2j+n-1)!/(2j)! / z^(2j+n) ]
static inline dd dd_psi_n_asym(const dd& z, int n) {
  const dd iz = 1.0 / z, iz2 = iz * iz;
  const dd izn = dd_powi(iz, n);
  dd s = GG2_FACT[n - 1] * izn + 0.5 * GG2_FACT[n] * (izn * iz);
  dd pw = izn;
  for (int j = 0; j < 17; j++) {
    pw = pw * iz2;
    const dd t = dd(GG2_BT[n - 1][j][0], GG2_BT[n - 1][j][1]) * pw;
    s = s + t;
    if (std::fabs(t.hi) < 1e-34 * std::fabs(s.hi)) break;
  }
  return (n % 2) ? s : -s;
}

// psi^(n)(x), n >= 1, by psi^(n)(x) = psi^(n)(x+1) + (-1)^(n+1) n! x^-(n+1)
static inline dd dd_psi_n(dd x, int n) {
  dd acc(0.0);
  while (x.hi < 40.0) {
    acc = acc + dd_powi(x, -(n + 1));
    x = x + 1.0;
  }
  const dd r = dd_psi_n_asym(x, n);
  const double sg = (n % 2) ? GG2_FACT[n] : -GG2_FACT[n];
  return r + sg * acc;
}

// psi(z) - log z, z >= 40: -1/(2z) - sum_j B_2j/(2j z^2j)
static inline dd dd_psi0_rest(const dd& z) {
  const dd iz = 1.0 / z, iz2 = iz * iz;
  dd s = -0.5 * iz, pw(1.0);
  for (int j = 0; j < 17; j++) {
    pw = pw * iz2;
    const dd t = dd(GG2_B0[j][0], GG2_B0[j][1]) * pw;
    s = s - t;
    if (std::fabs(t.hi) < 1e-34 * std::fabs(s.hi)) break;
  }
  return s;
}

// psi(k) - log k in double-double: psi(k) = psi(k+N) - sum_{i<N} 1/(k+i)
static inline dd dd_B0(dd k) {
  dd acc(0.0), x = k;
  while (x.hi < 40.0) {
    acc = acc + 1.0 / x;
    x = x + 1.0;
  }
  return dd_log(x / k) + dd_psi0_rest(x) - acc;
}

// R = e^z - sum_{i<8} z^i/i!, |z| < 1, by its series
static inline double gg2_R(double z) {
  double t = 1.0;
  for (int i = 1; i <= 8; i++) t *= z / i;
  double s = 0.0;
  for (int i = 8; i < 60; i++) {
    s += t;
    t *= z / (i + 1);
    if (std::fabs(t) < 1e-18 * std::fabs(s)) break;
  }
  return s;
}

static inline double pw(double x, int n) { return std::pow(x, n); }
static inline dd pw(const dd& x, int n) { return dd_powi(x, n); }
static inline double to_dbl(double x) { return x; }
static inline double to_dbl(const dd& x) { return x.hi; }

// below these, the coefficients are formed in double (measured,
// gengamma2_measure.R); above, in double-double
#define GG2_KC 10.0
#define GG2_DC 10.0

// psi^(n)(k1) - psi^(n)(k), both shifted by the same count so that the
// difference of the logarithms is log(z1/z) with z1 - z = k1 - k exactly
static inline dd dd_psi_diff(dd a, dd b, int n) {
  if (n >= 1) return dd_psi_n(b, n) - dd_psi_n(a, n);
  dd acc(0.0);
  while (a.hi < 40.0 || b.hi < 40.0) {
    acc = acc + (1.0 / a - 1.0 / b);
    a = a + 1.0;
    b = b + 1.0;
  }
  return dd_log(b / a) + dd_psi0_rest(b) - dd_psi0_rest(a) + acc;
}

// A0 = p (lgamma(k1) - lgamma(k)) - psi(k), so that X = p log(y/m) + A0.
// For k >= 20 through Stirling's series, where the logarithms of k and k1
// combine into log1p(1/d) and nothing large is differenced.
static inline double gg2_A0(double d, double p) {
  const double k = d / p, k1 = (d + 1.0) / p;
  if (k < 20.0) {
    if (d < 10.0) return p * (R::lgammafn(k1) - R::lgammafn(k)) - R::digamma(k);
    // h = 1/p small against k: A0 = sum_{j>=2} h^(j-1) psi^(j-1)(k)/j!,
    // whose terms fall by about 1/d
    const double h = 1.0 / p;
    double s = 0.0, hp = 1.0, fj = 1.0;
    for (int j = 2; j < 40; j++) {
      hp *= h; fj *= j;
      const double t = hp * R::psigamma(k, j - 1.0) / fj;
      s += t;
      if (std::fabs(t) < 1e-17 * std::fabs(s)) break;
    }
    return s;
  }
  // p GL = (d + 1 - p/2) log1p(1/d) + log k - 1 + p (S(k1) - S(k)),
  // psi(k) = log k + B(k); so A0 = (d + 1 - p/2) log1p(1/d) - 1
  //   + p (S(k1) - S(k)) - B(k)
  auto S = [](double z) {
    const double iz = 1.0 / z, iz2 = iz * iz;
    return iz * (1.0 / 12.0 + iz2 * (-1.0 / 360.0 + iz2 * (1.0 / 1260.0 +
           iz2 * (-1.0 / 1680.0 + iz2 * (1.0 / 1188.0)))));
  };
  auto B = [](double z) {
    const double iz = 1.0 / z, iz2 = iz * iz;
    return -0.5 * iz - iz2 * (1.0 / 12.0 + iz2 * (-1.0 / 120.0 + iz2 * (1.0 / 252.0 +
           iz2 * (-1.0 / 240.0 + iz2 * (1.0 / 132.0)))));
  };
  // (d + 1 - p/2) log1p(u) - 1 = [log1p(u)/u - 1] + (1 - p/2) log1p(u),
  // u = 1/d, with the bracket from its series where it cancels
  const double u = 1.0 / d, l = std::log1p(u);
  double s1;
  if (u < 0.1) {
    s1 = 0.0;
    double t = 1.0;
    for (int j = 1; j < 40; j++) {
      t *= -u;
      const double term = t / (j + 1.0);
      s1 += term;
      if (std::fabs(term) < 1e-18 * std::fabs(s1)) break;
    }
  } else {
    s1 = l / u - 1.0;
  }
  return s1 + (1.0 - 0.5 * p) * l + p * (S(k1) - S(k)) - B(k);
}

// B0 = psi(k) - log(k), so that Q = expm1(X + B0)
static inline double gg2_B0(double k) {
  if (k < 20.0) return R::digamma(k) - std::log(k);
  const double iz = 1.0 / k, iz2 = iz * iz;
  return -0.5 * iz - iz2 * (1.0 / 12.0 + iz2 * (-1.0 / 120.0 + iz2 * (1.0 / 252.0 +
         iz2 * (-1.0 / 240.0 + iz2 * (1.0 / 132.0 + iz2 * (-691.0 / 32760.0 +
         iz2 * (1.0 / 12.0)))))));
}

namespace d7 {

// The data quantities of one observation. With Z = pC + Z0, below |Z| = 1
// the components are polynomials in Z and R = e^Z - sum_{i<8} Z^i/i!, and
// above it in X = pC + A0 and Q = expm1(Z); the powers are formed once here
// and read by every component.
struct GG2Data { double Z, R, Q, X, x2; double z[10]; };

inline GG2Data gg2_data(double Z, double pC, double A0) {
  GG2Data D;
  D.Z = Z; D.R = 0.0; D.Q = 0.0; D.X = 0.0; D.x2 = 0.0;
  if (std::fabs(Z) < 1.0) {
    D.R = gg2_R(Z);
    for (int j = 2; j <= 9; ++j) D.z[j] = std::pow(Z, j);
  } else {
    D.X = pC + A0;
    D.Q = std::expm1(Z);
    D.x2 = std::pow(D.X, 2);
  }
  return D;
}

template <class T>
inline void gengamma2_coef_score_mean(const T& d, const T& p, const T* PAv, const T* Gv, const T& Bs, double* K) {
  (void) d; (void) p; (void) PAv; (void) Gv; (void) Bs;
  const T c0 = ((T(1.0) / T(5040.0)))*d;
  const T c1 = ((T(1.0) / T(720.0)))*d;
  const T c2 = ((T(1.0) / T(120.0)))*d;
  const T c3 = ((T(1.0) / T(24.0)))*d;
  const T c4 = ((T(1.0) / T(6.0)))*d;
  const T c5 = ((T(1.0) / T(2.0)))*d;
  K[0] = to_dbl(d);
  K[1] = to_dbl(c0);
  K[2] = to_dbl(c1);
  K[3] = to_dbl(c2);
  K[4] = to_dbl(c3);
  K[5] = to_dbl(c4);
  K[6] = to_dbl(c5);
  K[7] = to_dbl(d);
  K[8] = to_dbl(d);
}

template <class T>
inline void gengamma2_coef_score_d(const T& d, const T& p, const T* PAv, const T* Gv, const T& Bs, double* K) {
  (void) d; (void) p; (void) PAv; (void) Gv; (void) Bs;
  const T G0 = Gv[0];
  const T c0 = ((T(1.0) / T(5040.0)))*d;
  const T c1 = ((T(1.0) / T(720.0)))*d;
  const T c2 = ((T(1.0) / T(120.0)))*d;
  const T c3 = ((T(1.0) / T(24.0)))*d;
  const T c4 = ((T(1.0) / T(6.0)))*d;
  const T c5 = ((T(1.0) / T(2.0)))*d;
  const T c6 = pw(p, -1);
  const T c7 = G0*d;
  const T c8 = -c6*c7;
  const T c9 = G0*c6;
  const T c10 = c7 - 1.0;
  K[9] = to_dbl(c6);
  K[10] = to_dbl(c8);
  K[11] = to_dbl(-c0*c9);
  K[12] = to_dbl(-c1*c9);
  K[13] = to_dbl(-c2*c9);
  K[14] = to_dbl(-c3*c9);
  K[15] = to_dbl(-c4*c9);
  K[16] = to_dbl(-c5*c9);
  K[17] = to_dbl(-c10*c6);
  K[18] = to_dbl(c8);
  K[19] = to_dbl(-Bs*c6);
}

template <class T>
inline void gengamma2_coef_score_p(const T& d, const T& p, const T* PAv, const T* Gv, const T& Bs, double* K) {
  (void) d; (void) p; (void) PAv; (void) Gv; (void) Bs;
  const T G0 = Gv[0];
  const T c0 = ((T(1.0) / T(5040.0)))*d;
  const T c1 = ((T(1.0) / T(720.0)))*d;
  const T c2 = ((T(1.0) / T(120.0)))*d;
  const T c3 = ((T(1.0) / T(24.0)))*d;
  const T c4 = ((T(1.0) / T(6.0)))*d;
  const T c5 = ((T(1.0) / T(2.0)))*d;
  const T c6 = pw(p, -1);
  const T c7 = G0*d;
  const T c10 = c7 - 1.0;
  const T c11 = pw(p, -2);
  const T c12 = c11*d;
  const T c13 = -c12;
  const T c14 = c0*c11;
  const T c15 = Bs + G0;
  const T c16 = c15 + c7;
  K[20] = to_dbl(c13);
  K[21] = to_dbl(c13);
  K[22] = to_dbl(c11*c7*(d + 1.0));
  K[23] = to_dbl(c6);
  K[24] = to_dbl(-c14);
  K[25] = to_dbl(c14*(c16 - 7.0));
  K[26] = to_dbl(c1*c11*(c16 - 6.0));
  K[27] = to_dbl(c11*c2*(c16 - 5.0));
  K[28] = to_dbl(c11*c3*(c16 - 4.0));
  K[29] = to_dbl(c11*c4*(c16 - 3.0));
  K[30] = to_dbl(c11*c5*(c16 - 2.0));
  K[31] = to_dbl(c13);
  K[32] = to_dbl(c12*(c10 + c15));
  K[33] = to_dbl(c12*c16);
  K[34] = to_dbl(c11*(Bs*d + p));
}

template <class T>
inline void gengamma2_coef_hess_mean_mean(const T& d, const T& p, const T* PAv, const T* Gv, const T& Bs, double* K) {
  (void) d; (void) p; (void) PAv; (void) Gv; (void) Bs;
  const T c0 = d*(p + 1.0);
  const T c1 = -c0;
  const T c2 = d*p;
  const T c3 = -c2;
  K[0] = to_dbl(c1);
  K[1] = to_dbl(c3);
  K[2] = to_dbl(-(T(1.0) / T(5040.0))*c0);
  K[3] = to_dbl(-(T(1.0) / T(720.0))*c0);
  K[4] = to_dbl(-(T(1.0) / T(120.0))*c0);
  K[5] = to_dbl(-(T(1.0) / T(24.0))*c0);
  K[6] = to_dbl(-(T(1.0) / T(6.0))*c0);
  K[7] = to_dbl(-(T(1.0) / T(2.0))*c0);
  K[8] = to_dbl(c1);
  K[9] = to_dbl(c1);
  K[10] = to_dbl(c3);
}

template <class T>
inline void gengamma2_coef_hess_d_d(const T& d, const T& p, const T* PAv, const T* Gv, const T& Bs, double* K) {
  (void) d; (void) p; (void) PAv; (void) Gv; (void) Bs;
  const T PA1 = PAv[1];
  const T G0 = Gv[0];
  const T G1 = Gv[1];
  const T c2 = d*p;
  const T c4 = pw(p, 2);
  const T c5 = pw(c4, -1);
  const T c6 = pw(G0, 2);
  const T c7 = c6*p;
  const T c8 = G1 + c7;
  const T c9 = c5*c8*d;
  const T c10 = -c9;
  const T c11 = 2.0*G0;
  const T c12 = c11*p;
  const T c13 = -c12;
  const T c14 = c2*c6;
  const T c15 = PA1 + c14;
  const T c16 = c13 + c15;
  const T c17 = -c16*c5;
  K[11] = to_dbl(c10);
  K[12] = to_dbl(c17);
  K[13] = to_dbl(-(T(1.0) / T(5040.0))*c9);
  K[14] = to_dbl(-(T(1.0) / T(720.0))*c9);
  K[15] = to_dbl(-(T(1.0) / T(120.0))*c9);
  K[16] = to_dbl(-(T(1.0) / T(24.0))*c9);
  K[17] = to_dbl(-(T(1.0) / T(6.0))*c9);
  K[18] = to_dbl(-(T(1.0) / T(2.0))*c9);
  K[19] = to_dbl(c10);
  K[20] = to_dbl(c10);
  K[21] = to_dbl(c17);
}

template <class T>
inline void gengamma2_coef_hess_p_p(const T& d, const T& p, const T* PAv, const T* Gv, const T& Bs, double* K) {
  (void) d; (void) p; (void) PAv; (void) Gv; (void) Bs;
  const T PA1 = PAv[1];
  const T G0 = Gv[0];
  const T G1 = Gv[1];
  const T c2 = d*p;
  const T c4 = pw(p, 2);
  const T c6 = pw(G0, 2);
  const T c7 = c6*p;
  const T c8 = G1 + c7;
  const T c11 = 2.0*G0;
  const T c12 = c11*p;
  const T c13 = -c12;
  const T c14 = c2*c6;
  const T c18 = pw(p, -3);
  const T c19 = c18*d;
  const T c20 = -c19;
  const T c21 = c11*c19*(d + 1.0);
  const T c22 = pw(d, 2);
  const T c23 = c22*c7;
  const T c24 = G1*d;
  const T c25 = PA1*d;
  const T c26 = G1*c22 + PA1 + 2.0*c14 + c23 + 2.0*c24 + 2.0*c25 + c8;
  const T c27 = pw(p, -4);
  const T c28 = c27*d;
  const T c29 = -c11*c2;
  const T c30 = c14 + c29;
  const T c31 = PA1*c22 - c12*c22 + 2.0*c23 + c30 + c4 + c7*pw(d, 3);
  const T c32 = ((T(1.0) / T(5040.0)))*c19;
  const T c33 = 2.0*Bs;
  const T c34 = Bs*p;
  const T c35 = G0*p;
  const T c36 = G0*c2;
  const T c37 = c33*c36;
  const T c38 = pw(Bs, 2);
  const T c39 = c33*c35;
  const T c40 = c26 + c37 + c38*p + c39;
  const T c41 = 12.0*p;
  const T c42 = 6.0*p;
  const T c43 = -G0*c42;
  const T c44 = -4.0*c35;
  const T c45 = G0*d;
  const T c46 = Bs + G0;
  const T c47 = c45 + c46;
  K[22] = to_dbl(c20);
  K[23] = to_dbl(c20);
  K[24] = to_dbl(c21);
  K[25] = to_dbl(c21);
  K[26] = to_dbl(-c26*c28);
  K[27] = to_dbl(-c27*c31);
  K[28] = to_dbl(-c32);
  K[29] = to_dbl(c32*(c11*d + c11 + c33 - 7.0));
  K[30] = to_dbl(-(T(1.0) / T(5040.0))*c28*(-14.0*c34 - 14.0*c35 - 14.0*c36 + c40 + 42.0*p));
  K[31] = to_dbl(-(T(1.0) / T(720.0))*c28*(-Bs*c41 - G0*c41 - 12.0*c36 + c40 + 30.0*p));
  K[32] = to_dbl(-(T(1.0) / T(120.0))*c28*(-10.0*c34 - 10.0*c35 - 10.0*c36 + c40 + 20.0*p));
  K[33] = to_dbl(-(T(1.0) / T(24.0))*c28*(-8.0*c34 - 8.0*c35 - 8.0*c36 + c40 + c41));
  K[34] = to_dbl(-(T(1.0) / T(6.0))*c28*(-Bs*c42 - 6.0*c36 + c40 + c42 + c43));
  K[35] = to_dbl(c20);
  K[36] = to_dbl(-(T(1.0) / T(2.0))*c28*(-4.0*c34 - 4.0*c36 + c40 + c44 + 2.0*p));
  K[37] = to_dbl(2.0*c19*c47);
  K[38] = to_dbl(-c28*(c13 + c29 - c33*p + c40));
  K[39] = to_dbl(-c28*c40);
  K[40] = to_dbl(-c27*(c2*c38 + c22*c39 + c31 + c37));
}

template <class T>
inline void gengamma2_coef_hessian_rest(const T& d, const T& p, const T* PAv, const T* Gv, const T& Bs, double* K) {
  (void) d; (void) p; (void) PAv; (void) Gv; (void) Bs;
  const T PA1 = PAv[1];
  const T G0 = Gv[0];
  const T G1 = Gv[1];
  const T c2 = d*p;
  const T c4 = pw(p, 2);
  const T c5 = pw(c4, -1);
  const T c6 = pw(G0, 2);
  const T c7 = c6*p;
  const T c8 = G1 + c7;
  const T c11 = 2.0*G0;
  const T c12 = c11*p;
  const T c13 = -c12;
  const T c14 = c2*c6;
  const T c15 = PA1 + c14;
  const T c16 = c13 + c15;
  const T c18 = pw(p, -3);
  const T c19 = c18*d;
  const T c22 = pw(d, 2);
  const T c23 = c22*c7;
  const T c24 = G1*d;
  const T c25 = PA1*d;
  const T c29 = -c11*c2;
  const T c30 = c14 + c29;
  const T c32 = ((T(1.0) / T(5040.0)))*c19;
  const T c35 = G0*p;
  const T c36 = G0*c2;
  const T c42 = 6.0*p;
  const T c43 = -G0*c42;
  const T c44 = -4.0*c35;
  const T c45 = G0*d;
  const T c46 = Bs + G0;
  const T c47 = c45 + c46;
  const T c48 = c45 - 1.0;
  const T c49 = ((T(1.0) / T(5040.0)))*c45;
  const T c50 = d/p;
  const T c51 = -c50*(G0 + c48);
  const T c52 = ((T(1.0) / T(5040.0)))*c50;
  const T c53 = -c50*(c46 + c48);
  const T c54 = -c45*c5;
  const T c55 = c24 + c8;
  const T c56 = c15 + c55;
  const T c57 = -c35;
  const T c58 = c23 + c25 + c30 + c57;
  const T c59 = Bs*c35;
  const T c60 = c56 + c59;
  K[41] = to_dbl(c45);
  K[42] = to_dbl(c48);
  K[43] = to_dbl(c49);
  K[44] = to_dbl(((T(1.0) / T(720.0)))*c45);
  K[45] = to_dbl(((T(1.0) / T(120.0)))*c45);
  K[46] = to_dbl(((T(1.0) / T(24.0)))*c45);
  K[47] = to_dbl(((T(1.0) / T(6.0)))*c45);
  K[48] = to_dbl(((T(1.0) / T(2.0)))*c45);
  K[49] = to_dbl(c45);
  K[50] = to_dbl(c45);
  K[51] = to_dbl(c48);
  K[52] = to_dbl(c50);
  K[53] = to_dbl(c50);
  K[54] = to_dbl(c51);
  K[55] = to_dbl(c51);
  K[56] = to_dbl(c52);
  K[57] = to_dbl(-c52*(c47 - 8.0));
  K[58] = to_dbl(-(T(1.0) / T(720.0))*c50*(c47 - 7.0));
  K[59] = to_dbl(-(T(1.0) / T(120.0))*c50*(c47 - 6.0));
  K[60] = to_dbl(-(T(1.0) / T(24.0))*c50*(c47 - 5.0));
  K[61] = to_dbl(-(T(1.0) / T(6.0))*c50*(c47 - 4.0));
  K[62] = to_dbl(-(T(1.0) / T(2.0))*c50*(c47 - 3.0));
  K[63] = to_dbl(c50);
  K[64] = to_dbl(-c50*(c47 - 2.0));
  K[65] = to_dbl(c53);
  K[66] = to_dbl(c53);
  K[67] = to_dbl(c54);
  K[68] = to_dbl(c54);
  K[69] = to_dbl(c19*c56);
  K[70] = to_dbl(c18*c58);
  K[71] = to_dbl(-c49*c5);
  K[72] = to_dbl(c32*(-7.0*c35 + c60));
  K[73] = to_dbl(((T(1.0) / T(720.0)))*c19*(c43 + c60));
  K[74] = to_dbl(((T(1.0) / T(120.0)))*c19*(-5.0*c35 + c60));
  K[75] = to_dbl(((T(1.0) / T(24.0)))*c19*(c44 + c60));
  K[76] = to_dbl(((T(1.0) / T(6.0)))*c19*(-3.0*c35 + c60));
  K[77] = to_dbl(((T(1.0) / T(2.0)))*c19*(c16 + c55 + c59));
  K[78] = to_dbl(c54);
  K[79] = to_dbl(c19*(c57 + c60));
  K[80] = to_dbl(c19*c60);
  K[81] = to_dbl(c18*(Bs*c36 + c58));
}

template <class T>
inline void gengamma2_coef_expected_mean_mean(const T& d, const T& p, const T* PAv, const T* Gv, const T& Bs, double* K) {
  (void) d; (void) p; (void) PAv; (void) Gv; (void) Bs;
  const T c0 = d*p;
  K[0] = to_dbl(-c0);
}

template <class T>
inline void gengamma2_coef_expected_d_d(const T& d, const T& p, const T* PAv, const T* Gv, const T& Bs, double* K) {
  (void) d; (void) p; (void) PAv; (void) Gv; (void) Bs;
  const T PA1 = PAv[1];
  const T G0 = Gv[0];
  const T c0 = d*p;
  const T c1 = pw(p, 2);
  const T c2 = G0*p;
  const T c3 = pw(G0, 2);
  const T c4 = c0*c3;
  K[1] = to_dbl(-(PA1 - 2.0*c2 + c4)/c1);
}

template <class T>
inline void gengamma2_coef_expected_p_p(const T& d, const T& p, const T* PAv, const T* Gv, const T& Bs, double* K) {
  (void) d; (void) p; (void) PAv; (void) Gv; (void) Bs;
  const T PA1 = PAv[1];
  const T G0 = Gv[0];
  const T c0 = d*p;
  const T c1 = pw(p, 2);
  const T c3 = pw(G0, 2);
  const T c4 = c0*c3;
  const T c5 = pw(d, 2);
  const T c6 = G0*c1;
  const T c7 = c3*p;
  const T c8 = G0*d;
  const T c9 = G0*c5;
  const T c10 = c5*c7;
  const T c11 = -2.0*G0*c0 + c4;
  K[2] = to_dbl(-(PA1*c0 + PA1*c5 - 2.0*c1*c8 + c1 + 2.0*c10 + c11 - 2.0*c6 + c7*pw(d, 3) - 2.0*c9*p)/pw(p, 4));
}

template <class T>
inline void gengamma2_coef_expected_hessian_rest(const T& d, const T& p, const T* PAv, const T* Gv, const T& Bs, double* K) {
  (void) d; (void) p; (void) PAv; (void) Gv; (void) Bs;
  const T PA1 = PAv[1];
  const T G0 = Gv[0];
  const T c0 = d*p;
  const T c1 = pw(p, 2);
  const T c2 = G0*p;
  const T c3 = pw(G0, 2);
  const T c4 = c0*c3;
  const T c5 = pw(d, 2);
  const T c6 = G0*c1;
  const T c7 = c3*p;
  const T c8 = G0*d;
  const T c9 = G0*c5;
  const T c10 = c5*c7;
  const T c11 = -2.0*G0*c0 + c4;
  K[3] = to_dbl(c8 - 1.0);
  K[4] = to_dbl(-(c8 + c9 - d - p)/p);
  K[5] = to_dbl((PA1*d + c10 + c11 - c2 - c6)/pw(p, 3));
}

template <class T>
inline void gengamma2_coef_dexpected_mean_mean_mean(const T& d, const T& p, const T* PAv, const T* Gv, const T& Bs, double* K) {
  (void) d; (void) p; (void) PAv; (void) Gv; (void) Bs;
  const T c0 = 2.0*p;
  K[0] = to_dbl(c0*d);
}

template <class T>
inline void gengamma2_coef_dexpected_d_d_d(const T& d, const T& p, const T* PAv, const T* Gv, const T& Bs, double* K) {
  (void) d; (void) p; (void) PAv; (void) Gv; (void) Bs;
  const T PA2 = PAv[2];
  const T G0 = Gv[0];
  const T G1 = Gv[1];
  const T c3 = pw(p, 3);
  const T c4 = pw(c3, -1);
  const T c5 = pw(p, 2);
  const T c6 = pw(G0, 2)*c5;
  const T c7 = G1*p;
  const T c8 = 2.0*c7;
  const T c9 = G0*d;
  const T c10 = c8*c9;
  const T c11 = c10 - c8;
  K[4] = to_dbl(-c4*(PA2 + c11 + c6));
}

template <class T>
inline void gengamma2_coef_dexpected_p_p_p(const T& d, const T& p, const T* PAv, const T* Gv, const T& Bs, double* K) {
  (void) d; (void) p; (void) PAv; (void) Gv; (void) Bs;
  const T PA1 = PAv[1];
  const T PA2 = PAv[2];
  const T G0 = Gv[0];
  const T G1 = Gv[1];
  const T c0 = 2.0*p;
  const T c3 = pw(p, 3);
  const T c5 = pw(p, 2);
  const T c6 = pw(G0, 2)*c5;
  const T c7 = G1*p;
  const T c8 = 2.0*c7;
  const T c9 = G0*d;
  const T c10 = c8*c9;
  const T c13 = c6*d;
  const T c14 = PA1*p;
  const T c15 = 2.0*c14*c9;
  const T c17 = pw(d, 2);
  const T c18 = G1*c17;
  const T c19 = G0*p;
  const T c20 = 2.0*c19;
  const T c21 = 2.0*c5;
  const T c23 = G1*d;
  const T c24 = -c0*c23;
  const T c28 = PA1*c5;
  const T c29 = PA2*c17;
  const T c30 = 2.0*c3;
  const T c31 = G1*c5;
  const T c32 = -2.0*c31;
  const T c33 = c5*c9;
  const T c36 = PA1*d;
  const T c37 = c0*c36;
  const T c38 = 3.0*c6;
  const T c39 = pw(d, 3);
  const T c40 = G1*c20;
  const T c42 = 4.0*c18;
  const T c45 = G0*c3;
  const T c46 = c23*c5;
  const T c47 = 2.0*c17;
  const T c48 = 6.0*c19;
  const T c49 = PA1*c19;
  const T c50 = G0*c17;
  K[8] = to_dbl((G1*c39*c48 + PA1*c20*c39 + PA2*c39 + c10 + 3.0*c13 + c14*c47 + c15 + 4.0*c17*c49 + 6.0*c17*c6 - c18*c21 + c18*c48 + c24 + c28*d - 2.0*c28 + c29*p - 4.0*c3*c9 + c30 + c32 - 6.0*c33 - c37 + c38*c39 - c39*c8 + c40*pw(d, 4) - c42*p - 4.0*c45 - 4.0*c46 - 6.0*c5*c50)/pw(p, 6));
}

template <class T>
inline void gengamma2_coef_dexpected1_rest(const T& d, const T& p, const T* PAv, const T* Gv, const T& Bs, double* K) {
  (void) d; (void) p; (void) PAv; (void) Gv; (void) Bs;
  const T PA1 = PAv[1];
  const T PA2 = PAv[2];
  const T G0 = Gv[0];
  const T G1 = Gv[1];
  const T c0 = 2.0*p;
  const T c1 = -p;
  const T c2 = -d;
  const T c3 = pw(p, 3);
  const T c4 = pw(c3, -1);
  const T c5 = pw(p, 2);
  const T c6 = pw(G0, 2)*c5;
  const T c7 = G1*p;
  const T c8 = 2.0*c7;
  const T c9 = G0*d;
  const T c10 = c8*c9;
  const T c11 = c10 - c8;
  const T c12 = pw(p, -4);
  const T c13 = c6*d;
  const T c14 = PA1*p;
  const T c15 = 2.0*c14*c9;
  const T c16 = PA2*d;
  const T c17 = pw(d, 2);
  const T c18 = G1*c17;
  const T c19 = G0*p;
  const T c20 = 2.0*c19;
  const T c21 = 2.0*c5;
  const T c22 = G0*c21;
  const T c23 = G1*d;
  const T c24 = -c0*c23;
  const T c25 = -c22 + c24;
  const T c26 = c16 + c18*c20 + c25;
  const T c27 = pw(p, -5);
  const T c28 = PA1*c5;
  const T c29 = PA2*c17;
  const T c30 = 2.0*c3;
  const T c31 = G1*c5;
  const T c32 = -2.0*c31;
  const T c33 = c5*c9;
  const T c34 = 4.0*c33;
  const T c35 = c0*c18;
  const T c36 = PA1*d;
  const T c37 = c0*c36;
  const T c38 = 3.0*c6;
  const T c39 = pw(d, 3);
  const T c40 = G1*c20;
  const T c41 = c39*c40;
  const T c42 = 4.0*c18;
  const T c43 = c19*c42;
  const T c44 = c10 + c6;
  const T c45 = G0*c3;
  const T c46 = c23*c5;
  const T c47 = 2.0*c17;
  const T c49 = PA1*c19;
  const T c50 = G0*c17;
  const T c51 = pw(p, -1);
  const T c52 = c19 + c23;
  const T c53 = pw(c5, -1);
  const T c54 = G1 + PA1;
  const T c55 = c1 + c18;
  const T c56 = 2.0*c13;
  K[1] = to_dbl(c1);
  K[2] = to_dbl(c2);
  K[3] = to_dbl(0.0);
  K[5] = to_dbl(c12*(c11 + c13 + c15 + c26));
  K[6] = to_dbl(0.0);
  K[7] = to_dbl(-c27*(-G0*c30 + 4.0*c13 + c16*p + c17*c38 - c21*c23 + c25 + c28 + c29 + c32 - c34 - c35 + c37 + c41 + c43 + c44));
  K[9] = to_dbl(1.0 - c9);
  K[10] = to_dbl(c51*c52);
  K[11] = to_dbl(-c53*d*(c23 + c54));
  K[12] = to_dbl(c51*(c1 + c2 + c50 + c9));
  K[13] = to_dbl(-c53*(c0*c9 + c52 + c55));
  K[14] = to_dbl(c4*d*(c19 + 2.0*c23 + c36 + c54 + c55 + c9*p));
  K[15] = to_dbl(0.0);
  K[16] = to_dbl(c12*(c14 + c26 - c31 + c44 + c56 - c7));
  K[17] = to_dbl(c27*(-c10 + c14 - c15 + c22 + 3.0*c23*p + c28 - c29 + c31 + c34 + c35 - c36*p - c41 - c43 + c45 + c46 - c47*c49 - c47*c6 - c56 + c7));
}

inline void gengamma2_coefs_score_mean(double dv, double pv, double* K) {
  const double k0d = dv / pv;
  if (k0d < GG2_KC && dv < GG2_DC) {
    const double k1d = (dv + 1.0) / pv;
    double PAv[12] = {0}, Gv[12] = {0};
    (void) k1d;
    gengamma2_coef_score_mean<double>(dv, pv, PAv, Gv, 0.0, K);
  } else {
    const dd d(dv), p(pv);
    const dd k0 = d / p, k1 = (d + 1.0) / p;
    dd PAv[12], Gv[12];
    (void) k0; (void) k1;
    gengamma2_coef_score_mean<dd>(d, p, PAv, Gv, dd(0.0), K);
  }
}

inline void gengamma2_coefs_score_d(double dv, double pv, double* K) {
  const double k0d = dv / pv;
  if (k0d < GG2_KC && dv < GG2_DC) {
    const double k1d = (dv + 1.0) / pv;
    double PAv[12] = {0}, Gv[12] = {0};
    (void) k1d;
    Gv[0] = R::psigamma(k1d, 0.0) - R::psigamma(k0d, 0.0);
    gengamma2_coef_score_d<double>(dv, pv, PAv, Gv, gg2_B0(k0d), K);
  } else {
    const dd d(dv), p(pv);
    const dd k0 = d / p, k1 = (d + 1.0) / p;
    dd PAv[12], Gv[12];
    (void) k0; (void) k1;
    Gv[0] = dd_psi_diff(k0, k1, 0);
    gengamma2_coef_score_d<dd>(d, p, PAv, Gv, dd_B0(k0), K);
  }
}

inline void gengamma2_coefs_score_p(double dv, double pv, double* K) {
  const double k0d = dv / pv;
  if (k0d < GG2_KC && dv < GG2_DC) {
    const double k1d = (dv + 1.0) / pv;
    double PAv[12] = {0}, Gv[12] = {0};
    (void) k1d;
    Gv[0] = R::psigamma(k1d, 0.0) - R::psigamma(k0d, 0.0);
    gengamma2_coef_score_p<double>(dv, pv, PAv, Gv, gg2_B0(k0d), K);
  } else {
    const dd d(dv), p(pv);
    const dd k0 = d / p, k1 = (d + 1.0) / p;
    dd PAv[12], Gv[12];
    (void) k0; (void) k1;
    Gv[0] = dd_psi_diff(k0, k1, 0);
    gengamma2_coef_score_p<dd>(d, p, PAv, Gv, dd_B0(k0), K);
  }
}

inline void gengamma2_coefs_hess_mean_mean(double dv, double pv, double* K) {
  const double k0d = dv / pv;
  if (k0d < GG2_KC && dv < GG2_DC) {
    const double k1d = (dv + 1.0) / pv;
    double PAv[12] = {0}, Gv[12] = {0};
    (void) k1d;
    gengamma2_coef_hess_mean_mean<double>(dv, pv, PAv, Gv, 0.0, K);
  } else {
    const dd d(dv), p(pv);
    const dd k0 = d / p, k1 = (d + 1.0) / p;
    dd PAv[12], Gv[12];
    (void) k0; (void) k1;
    gengamma2_coef_hess_mean_mean<dd>(d, p, PAv, Gv, dd(0.0), K);
  }
}

inline void gengamma2_coefs_hess_d_d(double dv, double pv, double* K) {
  const double k0d = dv / pv;
  if (k0d < GG2_KC && dv < GG2_DC) {
    const double k1d = (dv + 1.0) / pv;
    double PAv[12] = {0}, Gv[12] = {0};
    (void) k1d;
    PAv[1] = R::psigamma(k0d, 1.0);
    Gv[0] = R::psigamma(k1d, 0.0) - R::psigamma(k0d, 0.0);
    Gv[1] = R::psigamma(k1d, 1.0) - R::psigamma(k0d, 1.0);
    gengamma2_coef_hess_d_d<double>(dv, pv, PAv, Gv, 0.0, K);
  } else {
    const dd d(dv), p(pv);
    const dd k0 = d / p, k1 = (d + 1.0) / p;
    dd PAv[12], Gv[12];
    (void) k0; (void) k1;
    PAv[1] = dd_psi_n(k0, 1);
    Gv[0] = dd_psi_diff(k0, k1, 0);
    Gv[1] = dd_psi_diff(k0, k1, 1);
    gengamma2_coef_hess_d_d<dd>(d, p, PAv, Gv, dd(0.0), K);
  }
}

inline void gengamma2_coefs_hess_p_p(double dv, double pv, double* K) {
  const double k0d = dv / pv;
  if (k0d < GG2_KC && dv < GG2_DC) {
    const double k1d = (dv + 1.0) / pv;
    double PAv[12] = {0}, Gv[12] = {0};
    (void) k1d;
    PAv[1] = R::psigamma(k0d, 1.0);
    Gv[0] = R::psigamma(k1d, 0.0) - R::psigamma(k0d, 0.0);
    Gv[1] = R::psigamma(k1d, 1.0) - R::psigamma(k0d, 1.0);
    gengamma2_coef_hess_p_p<double>(dv, pv, PAv, Gv, gg2_B0(k0d), K);
  } else {
    const dd d(dv), p(pv);
    const dd k0 = d / p, k1 = (d + 1.0) / p;
    dd PAv[12], Gv[12];
    (void) k0; (void) k1;
    PAv[1] = dd_psi_n(k0, 1);
    Gv[0] = dd_psi_diff(k0, k1, 0);
    Gv[1] = dd_psi_diff(k0, k1, 1);
    gengamma2_coef_hess_p_p<dd>(d, p, PAv, Gv, dd_B0(k0), K);
  }
}

inline void gengamma2_coefs_expected_mean_mean(double dv, double pv, double* K) {
  const double k0d = dv / pv;
  if (k0d < GG2_KC && dv < GG2_DC) {
    const double k1d = (dv + 1.0) / pv;
    double PAv[12] = {0}, Gv[12] = {0};
    (void) k1d;
    gengamma2_coef_expected_mean_mean<double>(dv, pv, PAv, Gv, 0.0, K);
  } else {
    const dd d(dv), p(pv);
    const dd k0 = d / p, k1 = (d + 1.0) / p;
    dd PAv[12], Gv[12];
    (void) k0; (void) k1;
    gengamma2_coef_expected_mean_mean<dd>(d, p, PAv, Gv, dd(0.0), K);
  }
}

inline void gengamma2_coefs_expected_d_d(double dv, double pv, double* K) {
  const double k0d = dv / pv;
  if (k0d < GG2_KC && dv < GG2_DC) {
    const double k1d = (dv + 1.0) / pv;
    double PAv[12] = {0}, Gv[12] = {0};
    (void) k1d;
    PAv[1] = R::psigamma(k0d, 1.0);
    Gv[0] = R::psigamma(k1d, 0.0) - R::psigamma(k0d, 0.0);
    gengamma2_coef_expected_d_d<double>(dv, pv, PAv, Gv, 0.0, K);
  } else {
    const dd d(dv), p(pv);
    const dd k0 = d / p, k1 = (d + 1.0) / p;
    dd PAv[12], Gv[12];
    (void) k0; (void) k1;
    PAv[1] = dd_psi_n(k0, 1);
    Gv[0] = dd_psi_diff(k0, k1, 0);
    gengamma2_coef_expected_d_d<dd>(d, p, PAv, Gv, dd(0.0), K);
  }
}

inline void gengamma2_coefs_expected_p_p(double dv, double pv, double* K) {
  const double k0d = dv / pv;
  if (k0d < GG2_KC && dv < GG2_DC) {
    const double k1d = (dv + 1.0) / pv;
    double PAv[12] = {0}, Gv[12] = {0};
    (void) k1d;
    PAv[1] = R::psigamma(k0d, 1.0);
    Gv[0] = R::psigamma(k1d, 0.0) - R::psigamma(k0d, 0.0);
    gengamma2_coef_expected_p_p<double>(dv, pv, PAv, Gv, 0.0, K);
  } else {
    const dd d(dv), p(pv);
    const dd k0 = d / p, k1 = (d + 1.0) / p;
    dd PAv[12], Gv[12];
    (void) k0; (void) k1;
    PAv[1] = dd_psi_n(k0, 1);
    Gv[0] = dd_psi_diff(k0, k1, 0);
    gengamma2_coef_expected_p_p<dd>(d, p, PAv, Gv, dd(0.0), K);
  }
}

inline void gengamma2_coefs_dexpected_mean_mean_mean(double dv, double pv, double* K) {
  const double k0d = dv / pv;
  if (k0d < GG2_KC && dv < GG2_DC) {
    const double k1d = (dv + 1.0) / pv;
    double PAv[12] = {0}, Gv[12] = {0};
    (void) k1d;
    gengamma2_coef_dexpected_mean_mean_mean<double>(dv, pv, PAv, Gv, 0.0, K);
  } else {
    const dd d(dv), p(pv);
    const dd k0 = d / p, k1 = (d + 1.0) / p;
    dd PAv[12], Gv[12];
    (void) k0; (void) k1;
    gengamma2_coef_dexpected_mean_mean_mean<dd>(d, p, PAv, Gv, dd(0.0), K);
  }
}

inline void gengamma2_coefs_dexpected_d_d_d(double dv, double pv, double* K) {
  const double k0d = dv / pv;
  if (k0d < GG2_KC && dv < GG2_DC) {
    const double k1d = (dv + 1.0) / pv;
    double PAv[12] = {0}, Gv[12] = {0};
    (void) k1d;
    PAv[2] = R::psigamma(k0d, 2.0);
    Gv[0] = R::psigamma(k1d, 0.0) - R::psigamma(k0d, 0.0);
    Gv[1] = R::psigamma(k1d, 1.0) - R::psigamma(k0d, 1.0);
    gengamma2_coef_dexpected_d_d_d<double>(dv, pv, PAv, Gv, 0.0, K);
  } else {
    const dd d(dv), p(pv);
    const dd k0 = d / p, k1 = (d + 1.0) / p;
    dd PAv[12], Gv[12];
    (void) k0; (void) k1;
    PAv[2] = dd_psi_n(k0, 2);
    Gv[0] = dd_psi_diff(k0, k1, 0);
    Gv[1] = dd_psi_diff(k0, k1, 1);
    gengamma2_coef_dexpected_d_d_d<dd>(d, p, PAv, Gv, dd(0.0), K);
  }
}

inline void gengamma2_coefs_dexpected_p_p_p(double dv, double pv, double* K) {
  const double k0d = dv / pv;
  if (k0d < GG2_KC && dv < GG2_DC) {
    const double k1d = (dv + 1.0) / pv;
    double PAv[12] = {0}, Gv[12] = {0};
    (void) k1d;
    PAv[1] = R::psigamma(k0d, 1.0);
    PAv[2] = R::psigamma(k0d, 2.0);
    Gv[0] = R::psigamma(k1d, 0.0) - R::psigamma(k0d, 0.0);
    Gv[1] = R::psigamma(k1d, 1.0) - R::psigamma(k0d, 1.0);
    gengamma2_coef_dexpected_p_p_p<double>(dv, pv, PAv, Gv, 0.0, K);
  } else {
    const dd d(dv), p(pv);
    const dd k0 = d / p, k1 = (d + 1.0) / p;
    dd PAv[12], Gv[12];
    (void) k0; (void) k1;
    PAv[1] = dd_psi_n(k0, 1);
    PAv[2] = dd_psi_n(k0, 2);
    Gv[0] = dd_psi_diff(k0, k1, 0);
    Gv[1] = dd_psi_diff(k0, k1, 1);
    gengamma2_coef_dexpected_p_p_p<dd>(d, p, PAv, Gv, dd(0.0), K);
  }
}

inline double gengamma2_score_mean(const double* K, const GG2Data& D, double im) {
  const double Z = D.Z, R = D.R, Q = D.Q, X = D.X;
  (void) Z; (void) R; (void) Q; (void) X; (void) im;
  if (std::fabs(Z) < 1.0) {
    const double u0 = D.z[7];
    const double u1 = D.z[6];
    const double u2 = D.z[5];
    const double u3 = D.z[4];
    const double u4 = D.z[3];
    const double u5 = D.z[2];
    return (K[1]*u0 + K[2]*u1 + K[3]*u2 + K[4]*u3 + K[5]*u4 + K[6]*u5 + K[7]*Z + K[8]*R) * im;
  }
  return (K[0]*Q) * im;
}

inline double gengamma2_score_d(const double* K, const GG2Data& D, double im) {
  const double Z = D.Z, R = D.R, Q = D.Q, X = D.X;
  (void) Z; (void) R; (void) Q; (void) X; (void) im;
  if (std::fabs(Z) < 1.0) {
    const double u0 = D.z[7];
    const double u1 = D.z[6];
    const double u2 = D.z[5];
    const double u3 = D.z[4];
    const double u4 = D.z[3];
    const double u5 = D.z[2];
    return (K[11]*u0 + K[12]*u1 + K[13]*u2 + K[14]*u3 + K[15]*u4 + K[16]*u5 + K[17]*Z + K[18]*R + K[19]);
  }
  return (K[10]*Q + K[9]*X);
}

inline double gengamma2_score_p(const double* K, const GG2Data& D, double im) {
  const double Z = D.Z, R = D.R, Q = D.Q, X = D.X;
  (void) Z; (void) R; (void) Q; (void) X; (void) im;
  if (std::fabs(Z) < 1.0) {
    const double u0 = D.z[7];
    const double u1 = D.z[6];
    const double u2 = D.z[5];
    const double u3 = D.z[4];
    const double u4 = D.z[3];
    const double u5 = D.z[2];
    return (K[24]*D.z[8] + K[25]*u0 + K[26]*u1 + K[27]*u2 + K[28]*u3 + K[29]*u4 + K[30]*u5 + K[31]*R*Z + K[32]*Z + K[33]*R + K[34]);
  }
  return (K[20]*Q*X + K[21]*X + K[22]*Q + K[23]);
}

inline double gengamma2_hess_mean_mean(const double* K, const GG2Data& D, double im) {
  const double Z = D.Z, R = D.R, Q = D.Q, X = D.X;
  (void) Z; (void) R; (void) Q; (void) X; (void) im;
  if (std::fabs(Z) < 1.0) {
    const double u0 = D.z[7];
    const double u1 = D.z[6];
    const double u2 = D.z[5];
    const double u3 = D.z[4];
    const double u4 = D.z[3];
    const double u5 = D.z[2];
    return (K[10] + K[2]*u0 + K[3]*u1 + K[4]*u2 + K[5]*u3 + K[6]*u4 + K[7]*u5 + K[8]*Z + K[9]*R) * std::pow(im, 2);
  }
  return (K[0]*Q + K[1]) * std::pow(im, 2);
}

inline double gengamma2_hess_d_d(const double* K, const GG2Data& D, double im) {
  const double Z = D.Z, R = D.R, Q = D.Q, X = D.X;
  (void) Z; (void) R; (void) Q; (void) X; (void) im;
  if (std::fabs(Z) < 1.0) {
    const double u0 = D.z[7];
    const double u1 = D.z[6];
    const double u2 = D.z[5];
    const double u3 = D.z[4];
    const double u4 = D.z[3];
    const double u5 = D.z[2];
    return (K[13]*u0 + K[14]*u1 + K[15]*u2 + K[16]*u3 + K[17]*u4 + K[18]*u5 + K[19]*Z + K[20]*R + K[21]);
  }
  return (K[11]*Q + K[12]);
}

inline double gengamma2_hess_p_p(const double* K, const GG2Data& D, double im) {
  const double Z = D.Z, R = D.R, Q = D.Q, X = D.X;
  (void) Z; (void) R; (void) Q; (void) X; (void) im;
  if (std::fabs(Z) < 1.0) {
    const double u0 = D.z[7];
    const double u1 = D.z[6];
    const double u2 = D.z[5];
    const double u3 = D.z[4];
    const double u4 = D.z[3];
    const double u5 = D.z[2];
    const double u6 = R*Z;
    const double u7 = D.z[8];
    return (K[28]*D.z[9] + K[29]*u7 + K[30]*u0 + K[31]*u1 + K[32]*u2 + K[33]*u3 + K[34]*u4 + K[35]*R*u5 + K[36]*u5 + K[37]*u6 + K[38]*Z + K[39]*R + K[40]);
  }
  const double u0 = Q*X;
  const double u1 = D.x2;
  return (K[22]*Q*u1 + K[23]*u1 + K[24]*u0 + K[25]*X + K[26]*Q + K[27]);
}

inline double gengamma2_expected_mean_mean(const double* K, double im) {
  return (K[0]) * std::pow(im, 2);
}

inline double gengamma2_expected_d_d(const double* K, double im) {
  return (K[1]);
}

inline double gengamma2_expected_p_p(const double* K, double im) {
  return (K[2]);
}

inline double gengamma2_dexpected_mean_mean_mean(const double* K, double im) {
  return (K[0]) * std::pow(im, 3);
}

inline double gengamma2_dexpected_d_d_d(const double* K, double im) {
  return (K[4]);
}

inline double gengamma2_dexpected_p_p_p(const double* K, double im) {
  return (K[8]);
}

inline void gengamma2_score_curv(int k, double y, const double* th,
                                 double* out) {
  const double mv = th[0], dv = th[1], pv = th[2];
  const double A0 = gg2_A0(dv, pv);
  const double Z0 = A0 + gg2_B0(dv / pv);
  const double im = 1.0 / mv;
  const double dy = y - mv;
  const double pC = pv * (std::fabs(dy) < 0.5 * mv ? std::log1p(dy / mv) : std::log(y / mv));
  const double Z = pC + Z0;
  const GG2Data D = gg2_data(Z, pC, A0);
  double Kg[35], Kh[82];
  if (k == 0) {
    gengamma2_coefs_score_mean(dv, pv, Kg);
    gengamma2_coefs_hess_mean_mean(dv, pv, Kh);
    out[0] = gengamma2_score_mean(Kg, D, im);
    out[1] = gengamma2_hess_mean_mean(Kh, D, im);
  } else if (k == 1) {
    gengamma2_coefs_score_d(dv, pv, Kg);
    gengamma2_coefs_hess_d_d(dv, pv, Kh);
    out[0] = gengamma2_score_d(Kg, D, im);
    out[1] = gengamma2_hess_d_d(Kh, D, im);
  } else {
    gengamma2_coefs_score_p(dv, pv, Kg);
    gengamma2_coefs_hess_p_p(dv, pv, Kh);
    out[0] = gengamma2_score_p(Kg, D, im);
    out[1] = gengamma2_hess_p_p(Kh, D, im);
  }
}

inline void gengamma2_info_dinfo(int k, double y, const double* th,
                                 double* out) {
  const double mv = th[0], dv = th[1], pv = th[2];
  const double im = 1.0 / mv;
  double Ke[6], Kd[18];
  if (k == 0) {
    gengamma2_coefs_expected_mean_mean(dv, pv, Ke);
    gengamma2_coefs_dexpected_mean_mean_mean(dv, pv, Kd);
    out[0] = gengamma2_expected_mean_mean(Ke, im);
    out[1] = gengamma2_dexpected_mean_mean_mean(Kd, im);
  } else if (k == 1) {
    gengamma2_coefs_expected_d_d(dv, pv, Ke);
    gengamma2_coefs_dexpected_d_d_d(dv, pv, Kd);
    out[0] = gengamma2_expected_d_d(Ke, im);
    out[1] = gengamma2_dexpected_d_d_d(Kd, im);
  } else {
    gengamma2_coefs_expected_p_p(dv, pv, Ke);
    gengamma2_coefs_dexpected_p_p_p(dv, pv, Kd);
    out[0] = gengamma2_expected_p_p(Ke, im);
    out[1] = gengamma2_dexpected_p_p_p(Kd, im);
  }
}

// --- The joint registry ----------------------------------------------------
// The coefficient arrays of the kernels, set up as the kernels set them up:
// the gradient, the Hessian, the expected Hessian and its first derivatives
// each take the polygammas of (k0, k1) and the same Bs as in their kernel,
// in double below GG2_KC and GG2_DC and in double-double above them. The
// arrays are the K of the kernel, and the functions below read them.
inline void gengamma2_joint_coefs_grad(double dv, double pv, double* K) {
  const double k0d = dv / pv;
  if (k0d < GG2_KC && dv < GG2_DC) {
    const double k1d = (dv + 1.0) / pv;
    double PAv[12] = {0}, Gv[12] = {0};
    Gv[0] = R::psigamma(k1d, 0.0) - R::psigamma(k0d, 0.0);
    gengamma2_coef_score_mean<double>(dv, pv, PAv, Gv, gg2_B0(k0d), K);
    gengamma2_coef_score_d<double>(dv, pv, PAv, Gv, gg2_B0(k0d), K);
    gengamma2_coef_score_p<double>(dv, pv, PAv, Gv, gg2_B0(k0d), K);
  } else {
    const dd d(dv), p(pv);
    const dd k0 = d / p, k1 = (d + 1.0) / p;
    dd PAv[12], Gv[12];
    Gv[0] = dd_psi_diff(k0, k1, 0);
    gengamma2_coef_score_mean<dd>(d, p, PAv, Gv, dd_B0(k0), K);
    gengamma2_coef_score_d<dd>(d, p, PAv, Gv, dd_B0(k0), K);
    gengamma2_coef_score_p<dd>(d, p, PAv, Gv, dd_B0(k0), K);
  }
}

inline void gengamma2_joint_coefs_hess(double dv, double pv, double* K) {
  const double k0d = dv / pv;
  if (k0d < GG2_KC && dv < GG2_DC) {
    const double k1d = (dv + 1.0) / pv;
    double PAv[12] = {0}, Gv[12] = {0};
    PAv[1] = R::psigamma(k0d, 1.0);
    Gv[0] = R::psigamma(k1d, 0.0) - R::psigamma(k0d, 0.0);
    Gv[1] = R::psigamma(k1d, 1.0) - R::psigamma(k0d, 1.0);
    gengamma2_coef_hess_mean_mean<double>(dv, pv, PAv, Gv, gg2_B0(k0d), K);
    gengamma2_coef_hess_d_d<double>(dv, pv, PAv, Gv, gg2_B0(k0d), K);
    gengamma2_coef_hess_p_p<double>(dv, pv, PAv, Gv, gg2_B0(k0d), K);
    gengamma2_coef_hessian_rest<double>(dv, pv, PAv, Gv, gg2_B0(k0d), K);
  } else {
    const dd d(dv), p(pv);
    const dd k0 = d / p, k1 = (d + 1.0) / p;
    dd PAv[12], Gv[12];
    PAv[1] = dd_psi_n(k0, 1);
    Gv[0] = dd_psi_diff(k0, k1, 0);
    Gv[1] = dd_psi_diff(k0, k1, 1);
    gengamma2_coef_hess_mean_mean<dd>(d, p, PAv, Gv, dd_B0(k0), K);
    gengamma2_coef_hess_d_d<dd>(d, p, PAv, Gv, dd_B0(k0), K);
    gengamma2_coef_hess_p_p<dd>(d, p, PAv, Gv, dd_B0(k0), K);
    gengamma2_coef_hessian_rest<dd>(d, p, PAv, Gv, dd_B0(k0), K);
  }
}

inline void gengamma2_joint_coefs_expected(double dv, double pv, double* K) {
  const double k0d = dv / pv;
  if (k0d < GG2_KC && dv < GG2_DC) {
    const double k1d = (dv + 1.0) / pv;
    double PAv[12] = {0}, Gv[12] = {0};
    PAv[1] = R::psigamma(k0d, 1.0);
    Gv[0] = R::psigamma(k1d, 0.0) - R::psigamma(k0d, 0.0);
    gengamma2_coef_expected_mean_mean<double>(dv, pv, PAv, Gv, 0.0, K);
    gengamma2_coef_expected_d_d<double>(dv, pv, PAv, Gv, 0.0, K);
    gengamma2_coef_expected_p_p<double>(dv, pv, PAv, Gv, 0.0, K);
    gengamma2_coef_expected_hessian_rest<double>(dv, pv, PAv, Gv, 0.0, K);
  } else {
    const dd d(dv), p(pv);
    const dd k0 = d / p, k1 = (d + 1.0) / p;
    dd PAv[12], Gv[12];
    PAv[1] = dd_psi_n(k0, 1);
    Gv[0] = dd_psi_diff(k0, k1, 0);
    gengamma2_coef_expected_mean_mean<dd>(d, p, PAv, Gv, dd(0.0), K);
    gengamma2_coef_expected_d_d<dd>(d, p, PAv, Gv, dd(0.0), K);
    gengamma2_coef_expected_p_p<dd>(d, p, PAv, Gv, dd(0.0), K);
    gengamma2_coef_expected_hessian_rest<dd>(d, p, PAv, Gv, dd(0.0), K);
  }
}

inline void gengamma2_joint_coefs_dexpected(double dv, double pv, double* K) {
  const double k0d = dv / pv;
  if (k0d < GG2_KC && dv < GG2_DC) {
    const double k1d = (dv + 1.0) / pv;
    double PAv[12] = {0}, Gv[12] = {0};
    PAv[1] = R::psigamma(k0d, 1.0);
    PAv[2] = R::psigamma(k0d, 2.0);
    Gv[0] = R::psigamma(k1d, 0.0) - R::psigamma(k0d, 0.0);
    Gv[1] = R::psigamma(k1d, 1.0) - R::psigamma(k0d, 1.0);
    gengamma2_coef_dexpected_mean_mean_mean<double>(dv, pv, PAv, Gv, 0.0, K);
    gengamma2_coef_dexpected_d_d_d<double>(dv, pv, PAv, Gv, 0.0, K);
    gengamma2_coef_dexpected_p_p_p<double>(dv, pv, PAv, Gv, 0.0, K);
    gengamma2_coef_dexpected1_rest<double>(dv, pv, PAv, Gv, 0.0, K);
  } else {
    const dd d(dv), p(pv);
    const dd k0 = d / p, k1 = (d + 1.0) / p;
    dd PAv[12], Gv[12];
    PAv[1] = dd_psi_n(k0, 1);
    PAv[2] = dd_psi_n(k0, 2);
    Gv[0] = dd_psi_diff(k0, k1, 0);
    Gv[1] = dd_psi_diff(k0, k1, 1);
    gengamma2_coef_dexpected_mean_mean_mean<dd>(d, p, PAv, Gv, dd(0.0), K);
    gengamma2_coef_dexpected_d_d_d<dd>(d, p, PAv, Gv, dd(0.0), K);
    gengamma2_coef_dexpected_p_p_p<dd>(d, p, PAv, Gv, dd(0.0), K);
    gengamma2_coef_dexpected1_rest<dd>(d, p, PAv, Gv, dd(0.0), K);
  }
}

// The off-diagonal Hessian entries, from gengamma2_hessian_cpp: the branch
// |Z| < 1 polynomial in Z and R, and the branch above it in X and Q.
inline double gengamma2_hess_mean_d(const double* K, const GG2Data& D, double im) {
  const double Z = D.Z, R = D.R, Q = D.Q, X = D.X;
  (void) Z; (void) R; (void) Q; (void) X;
  if (std::fabs(Z) < 1.0) {
    const double u0 = D.z[7], u1 = D.z[6], u2 = D.z[5], u3 = D.z[4], u4 = D.z[3], u5 = D.z[2];
    return (K[43]*u0 + K[44]*u1 + K[45]*u2 + K[46]*u3 + K[47]*u4 + K[48]*u5 + K[49]*Z + K[50]*R + K[51]) * im;
  }
  return (K[41]*Q + K[42]) * im;
}

inline double gengamma2_hess_mean_p(const double* K, const GG2Data& D, double im) {
  const double Z = D.Z, R = D.R, Q = D.Q, X = D.X;
  (void) Z; (void) R; (void) Q; (void) X;
  if (std::fabs(Z) < 1.0) {
    const double u0 = D.z[7], u1 = D.z[6], u2 = D.z[5], u3 = D.z[4], u4 = D.z[3], u5 = D.z[2];
    const double u6 = R*Z;
    const double u7 = D.z[8];
    return (K[56]*u7 + K[57]*u0 + K[58]*u1 + K[59]*u2 + K[60]*u3 + K[61]*u4 + K[62]*u5 + K[63]*u6 + K[64]*Z + K[65]*R + K[66]) * im;
  }
  const double u0 = Q*X;
  return (K[52]*u0 + K[53]*X + K[54]*Q + K[55]) * im;
}

inline double gengamma2_hess_d_p(const double* K, const GG2Data& D, double im) {
  const double Z = D.Z, R = D.R, Q = D.Q, X = D.X;
  (void) Z; (void) R; (void) Q; (void) X; (void) im;
  if (std::fabs(Z) < 1.0) {
    const double u0 = D.z[7], u1 = D.z[6], u2 = D.z[5], u3 = D.z[4], u4 = D.z[3], u5 = D.z[2];
    const double u6 = R*Z;
    const double u7 = D.z[8];
    return (K[71]*u7 + K[72]*u0 + K[73]*u1 + K[74]*u2 + K[75]*u3 + K[76]*u4 + K[77]*u5 + K[78]*u6 + K[79]*Z + K[80]*R + K[81]);
  }
  const double u0 = Q*X;
  return (K[67]*u0 + K[68]*X + K[69]*Q + K[70]);
}

// The expected information, from gengamma2_expected_hessian_cpp.
inline double gengamma2_expected_mean_d(const double* K, double im) {
  return (K[3]) * im;
}

inline double gengamma2_expected_mean_p(const double* K, double im) {
  return (K[4]) * im;
}

inline double gengamma2_expected_d_p(const double* K) {
  return (K[5]);
}

// The derivatives of the expected information, from gengamma2_dexpected1_cpp:
// the pair comes first and the differentiating parameter last.
inline double gengamma2_dexpected_mean_mean_d(const double* K, double im) {
  return (K[1]) * std::pow(im, 2);
}

inline double gengamma2_dexpected_mean_mean_p(const double* K, double im) {
  return (K[2]) * std::pow(im, 2);
}

inline double gengamma2_dexpected_d_d_mean(const double* K, double im) {
  return (K[3]) * im;
}

inline double gengamma2_dexpected_d_d_p(const double* K) {
  return (K[5]);
}

inline double gengamma2_dexpected_p_p_mean(const double* K, double im) {
  return (K[6]) * im;
}

inline double gengamma2_dexpected_p_p_d(const double* K) {
  return (K[7]);
}

inline double gengamma2_dexpected_mean_d_mean(const double* K, double im) {
  return (K[9]) * std::pow(im, 2);
}

inline double gengamma2_dexpected_mean_d_d(const double* K, double im) {
  return (K[10]) * im;
}

inline double gengamma2_dexpected_mean_d_p(const double* K, double im) {
  return (K[11]) * im;
}

inline double gengamma2_dexpected_mean_p_mean(const double* K, double im) {
  return (K[12]) * std::pow(im, 2);
}

inline double gengamma2_dexpected_mean_p_d(const double* K, double im) {
  return (K[13]) * im;
}

inline double gengamma2_dexpected_mean_p_p(const double* K, double im) {
  return (K[14]) * im;
}

inline double gengamma2_dexpected_d_p_mean(const double* K, double im) {
  return (K[15]) * im;
}

inline double gengamma2_dexpected_d_p_d(const double* K) {
  return (K[16]);
}

inline double gengamma2_dexpected_d_p_p(const double* K) {
  return (K[17]);
}

// The data of one observation for the joint routers, as gengamma2_hessian_cpp
// forms them: im = 1/mean, and the data quantities at Z = pC + Z0.
inline GG2Data gengamma2_joint_data(double y, double mv, double dv, double pv, double& im) {
  const double A0 = gg2_A0(dv, pv);
  const double Z0 = A0 + gg2_B0(dv / pv);
  im = 1.0 / mv;
  const double dy = y - mv;
  const double pC = pv * (std::fabs(dy) < 0.5 * mv ? std::log1p(dy / mv) : std::log(y / mv));
  const double Z = pC + Z0;
  return gg2_data(Z, pC, A0);
}

// The routers of the joint registry, one per order: every component at one
// observation, on the parameter scale, in (mean, d, p). The Hessian and the
// expected information fill a full 3 x 3 matrix and the derivative of the
// expected information a 3 x 3 x 3 array, column-major, out[a + 3 b + 9 c]
// being the derivative of E[l_ab] in parameter c.
inline void gengamma2_joint_grad(double y, const double* th, double* out) {
  const double mv = th[0], dv = th[1], pv = th[2];
  double im = 0.0;
  const GG2Data D = gengamma2_joint_data(y, mv, dv, pv, im);
  double K[35] = {0};
  gengamma2_joint_coefs_grad(dv, pv, K);
  out[0] = gengamma2_score_mean(K, D, im);
  out[1] = gengamma2_score_d(K, D, im);
  out[2] = gengamma2_score_p(K, D, im);
}

inline void gengamma2_joint_hess(double y, const double* th, double* out) {
  const double mv = th[0], dv = th[1], pv = th[2];
  double im = 0.0;
  const GG2Data D = gengamma2_joint_data(y, mv, dv, pv, im);
  double K[82] = {0};
  gengamma2_joint_coefs_hess(dv, pv, K);
  out[0] = gengamma2_hess_mean_mean(K, D, im);
  out[4] = gengamma2_hess_d_d(K, D, im);
  out[8] = gengamma2_hess_p_p(K, D, im);
  out[1] = out[3] = gengamma2_hess_mean_d(K, D, im);
  out[2] = out[6] = gengamma2_hess_mean_p(K, D, im);
  out[5] = out[7] = gengamma2_hess_d_p(K, D, im);
}

inline void gengamma2_joint_info(double y, const double* th, double* out) {
  const double mv = th[0], dv = th[1], pv = th[2];
  const double im = 1.0 / mv;
  double K[6] = {0};
  gengamma2_joint_coefs_expected(dv, pv, K);
  out[0] = gengamma2_expected_mean_mean(K, im);
  out[4] = gengamma2_expected_d_d(K, im);
  out[8] = gengamma2_expected_p_p(K, im);
  out[1] = out[3] = gengamma2_expected_mean_d(K, im);
  out[2] = out[6] = gengamma2_expected_mean_p(K, im);
  out[5] = out[7] = gengamma2_expected_d_p(K);
}

inline void gengamma2_joint_dinfo(double y, const double* th, double* out) {
  const double mv = th[0], dv = th[1], pv = th[2];
  const double im = 1.0 / mv;
  double K[18] = {0};
  gengamma2_joint_coefs_dexpected(dv, pv, K);
  for (int j = 0; j < 27; ++j) out[j] = 0.0;
  // c = 0 (mean)
  out[0] = gengamma2_dexpected_mean_mean_mean(K, im);
  out[4] = gengamma2_dexpected_d_d_mean(K, im);
  out[8] = gengamma2_dexpected_p_p_mean(K, im);
  out[1] = out[3] = gengamma2_dexpected_mean_d_mean(K, im);
  out[2] = out[6] = gengamma2_dexpected_mean_p_mean(K, im);
  out[5] = out[7] = gengamma2_dexpected_d_p_mean(K, im);
  // c = 1 (d)
  out[9] = gengamma2_dexpected_mean_mean_d(K, im);
  out[13] = gengamma2_dexpected_d_d_d(K, im);
  out[17] = gengamma2_dexpected_p_p_d(K);
  out[10] = out[12] = gengamma2_dexpected_mean_d_d(K, im);
  out[11] = out[15] = gengamma2_dexpected_mean_p_d(K, im);
  out[14] = out[16] = gengamma2_dexpected_d_p_d(K);
  // c = 2 (p)
  out[18] = gengamma2_dexpected_mean_mean_p(K, im);
  out[22] = gengamma2_dexpected_d_d_p(K);
  out[26] = gengamma2_dexpected_p_p_p(K, im);
  out[19] = out[21] = gengamma2_dexpected_mean_d_p(K, im);
  out[20] = out[24] = gengamma2_dexpected_mean_p_p(K, im);
  out[23] = out[25] = gengamma2_dexpected_d_p_p(K);
}

// --- The parts route (pt_parts.h) -------------------------------------------
// The coefficient generators of the third and fourth derivatives, copied
// from coef_gengamma2_deriv3_cpp and coef_gengamma2_deriv4_cpp of
// gengamma2.cpp (same templates, same operations), and the per-observation
// expressions of the kernels gengamma2_deriv3_cpp and gengamma2_deriv4_cpp.

template <class T>
inline void gengamma2_coef_deriv3(const T& d, const T& p, const T* PAv, const T* Gv, const T& Bs, double* K) {
  (void) d; (void) p; (void) PAv; (void) Gv; (void) Bs;
  const T PA1 = PAv[1];
  const T PA2 = PAv[2];
  const T G0 = Gv[0];
  const T G1 = Gv[1];
  const T G2 = Gv[2];
  const T c0 = 3.0*p;
  const T c1 = pw(p, 2);
  const T c2 = d*(c0 + c1 + 2.0);
  const T c3 = d*p*(p + 3.0);
  const T c4 = p + 1.0;
  const T c5 = G0*d;
  const T c6 = c4*c5;
  const T c7 = -c6;
  const T c8 = c5*p;
  const T c9 = c5 - 1.0;
  const T c10 = c8 + c9;
  const T c11 = -c10;
  const T c12 = pw(p, -1);
  const T c13 = c12*d;
  const T c14 = c13*c4;
  const T c15 = -c14;
  const T c16 = 2.0*p;
  const T c17 = G0*p;
  const T c18 = c10 - c16 + c17;
  const T c19 = c13*(G0 + c18);
  const T c20 = 9.0*p;
  const T c21 = Bs*p;
  const T c22 = Bs + G0;
  const T c23 = c22 + c5;
  const T c24 = c17 + c21 + c23 + c8;
  const T c25 = 8.0*p;
  const T c26 = 7.0*p;
  const T c27 = 6.0*p;
  const T c28 = 5.0*p;
  const T c29 = 4.0*p;
  const T c30 = c13*(c18 + c21 + c22);
  const T c31 = pw(G0, 2);
  const T c32 = c31*p;
  const T c33 = G1 + c32;
  const T c34 = c13*c33;
  const T c35 = c12*c5;
  const T c36 = G1*d;
  const T c37 = PA1 + c33;
  const T c38 = c32*d + c36 + c37;
  const T c39 = -c17 + c38;
  const T c40 = d/c1;
  const T c41 = -c39*c40;
  const T c42 = -8.0*c17;
  const T c43 = Bs*c17;
  const T c44 = c38 + c43;
  const T c45 = ((T(1.0) / T(5040.0)))*c40;
  const T c46 = -7.0*c17;
  const T c47 = 6.0*c17;
  const T c48 = -c47;
  const T c49 = 4.0*c17;
  const T c50 = -c49;
  const T c51 = G0*c0;
  const T c52 = G0*c16;
  const T c53 = -c52;
  const T c54 = -c40*(c39 + c43);
  const T c55 = G0 + c9;
  const T c56 = 2.0*c40;
  const T c57 = -c55*c56;
  const T c58 = c16*c31;
  const T c59 = 2.0*c36 + c58*d;
  const T c60 = pw(d, 2);
  const T c61 = G1*c60;
  const T c62 = 2.0*PA1;
  const T c63 = c32*c60 + c37 + c61 + c62*d;
  const T c64 = c59 + c63;
  const T c65 = -c16*c5 + c53;
  const T c66 = pw(p, 3);
  const T c67 = pw(c66, -1);
  const T c68 = c67*d;
  const T c69 = c68*(c64 + c65);
  const T c70 = 2.0*Bs;
  const T c71 = 2.0*G0;
  const T c72 = Bs*c52 + c59;
  const T c73 = pw(Bs, 2);
  const T c74 = Bs*c16;
  const T c75 = c5*c74 + c63 + c73*p;
  const T c76 = c72 + c75;
  const T c77 = 14.0*p;
  const T c78 = 30.0*p;
  const T c79 = 12.0*p;
  const T c80 = 6.0*c21;
  const T c81 = c68*(c65 - c74 + c76);
  const T c82 = G1*c51;
  const T c83 = pw(G0, 3)*c1;
  const T c84 = G2 + c82 + c83;
  const T c85 = c68*c84;
  const T c86 = -c85;
  const T c87 = c83*d;
  const T c88 = G1*c0;
  const T c89 = c5*c88 + c87;
  const T c90 = PA2 - c88;
  const T c91 = c89 + c90;
  const T c92 = -c67*c91;
  const T c93 = c33*c68;
  const T c94 = -c93;
  const T c95 = G1*p;
  const T c96 = PA2 + c84;
  const T c97 = G2*d;
  const T c98 = PA1*c52 + c97;
  const T c99 = c89 + c96 + c98;
  const T c100 = c95 + c99;
  const T c101 = pw(p, -4);
  const T c102 = c101*d;
  const T c103 = PA2*d;
  const T c104 = c0*c36;
  const T c105 = -c104;
  const T c106 = PA1*c16;
  const T c107 = c1*c71;
  const T c108 = -c107 + c89;
  const T c109 = G1*c16;
  const T c110 = -c109;
  const T c111 = c60*c83;
  const T c112 = c111 + c51*c61;
  const T c113 = c110 + c112;
  const T c114 = c103 + c105 + c106*c5 + c108 + c113;
  const T c115 = c1*c31;
  const T c116 = G1*c27;
  const T c117 = -c116;
  const T c118 = Bs*c115;
  const T c119 = G1*c21 + c118;
  const T c120 = c119 + c99;
  const T c121 = ((T(1.0) / T(5040.0)))*c102;
  const T c122 = 6.0*c115;
  const T c123 = -c122;
  const T c124 = -G1*c29;
  const T c125 = 4.0*c115;
  const T c126 = -c125;
  const T c127 = 3.0*c115;
  const T c128 = -c127;
  const T c129 = -c95;
  const T c130 = 2.0*c115;
  const T c131 = -c130;
  const T c132 = c118*d;
  const T c133 = c5*c67;
  const T c134 = -c133;
  const T c135 = 2.0*c102;
  const T c136 = c135*c38;
  const T c137 = c16*c36;
  const T c138 = G2*c60;
  const T c139 = PA1*c29;
  const T c140 = PA1*c51;
  const T c141 = c140 + c96;
  const T c142 = 2.0*c103 + c116*c5 + c138 + c139*c5 + c141 + 2.0*c87 + 2.0*c97;
  const T c143 = c112 + c142;
  const T c144 = c106 + c109 + c137 + c143;
  const T c145 = pw(p, -5);
  const T c146 = c145*d;
  const T c147 = PA2*c60;
  const T c148 = PA1*p;
  const T c149 = c1*c5;
  const T c150 = -c29*c36;
  const T c151 = PA1*c0;
  const T c152 = c151*c5;
  const T c153 = PA1*c60;
  const T c154 = c0*c61;
  const T c155 = pw(d, 3);
  const T c156 = c155*c83;
  const T c157 = c155*c82 + c156;
  const T c158 = -c154 + c157;
  const T c159 = c108 + 2.0*c111 + c129 + c147 - c148 - 4.0*c149 + c150 + c152 + c153*c49 + c158 + c47*c61;
  const T c160 = c46 + c72;
  const T c161 = 14.0*c115;
  const T c162 = G0*c1;
  const T c163 = 42.0*c162;
  const T c164 = Bs*c162;
  const T c165 = -c36*c79;
  const T c166 = PA1*c74;
  const T c167 = c115*c70;
  const T c168 = c167*d + c36*c74;
  const T c169 = G1*c74 + c162*c73 + c166 + c167 + c168;
  const T c170 = c143 + c169;
  const T c171 = PA1*c79;
  const T c172 = -G1*c79 - c171;
  const T c173 = 12.0*c115;
  const T c174 = -c173;
  const T c175 = 30.0*c162;
  const T c176 = c36*p;
  const T c177 = 12.0*c162;
  const T c178 = -Bs*c177 - c173*d;
  const T c179 = 10.0*c115;
  const T c180 = 8.0*c115;
  const T c181 = c27*c36;
  const T c182 = -c181;
  const T c183 = PA1*c27;
  const T c184 = c117 - c183;
  const T c185 = 6.0*c1;
  const T c186 = G0*c185;
  const T c187 = -Bs*c186 - c122*d;
  const T c188 = c123 + c186;
  const T c189 = c1*c73;
  const T c190 = -c102;
  const T c191 = 3.0*c5;
  const T c192 = c101*c191*(d + 1.0);
  const T c193 = 3.0*c146;
  const T c194 = -c193*c64;
  const T c195 = c183*d;
  const T c196 = G2*c155;
  const T c197 = 3.0*c97;
  const T c198 = 3.0*c138;
  const T c199 = 3.0*c103;
  const T c200 = 3.0*c147;
  const T c201 = 3.0*c87;
  const T c202 = c153*c47;
  const T c203 = G1*c20;
  const T c204 = c203*c5;
  const T c205 = PA1*c20;
  const T c206 = c205*c5;
  const T c207 = 9.0*c17;
  const T c208 = 3.0*c111 + c207*c61;
  const T c209 = c196 + c197 + c198 + c199 + c200 + c201 + c202 + c204 + c206 + c208;
  const T c210 = c141 + c157;
  const T c211 = c209 + c210;
  const T c212 = c151 + c154 + c181 + c195 + c211 + c88;
  const T c213 = pw(p, -6);
  const T c214 = c213*d;
  const T c215 = pw(d, 4);
  const T c216 = c185*c5;
  const T c217 = c155*c88;
  const T c218 = c151*d;
  const T c219 = c208 - c27*c61;
  const T c220 = G1*c155*c207 + PA1*c155*c47 + PA2*c155 + c105 + c152 + c153*c207 + 3.0*c156 - c186*c60 + c215*c82 + c215*c83 - c216 - c217 - c218 + c219 + 2.0*c66 + c89;
  const T c221 = 3.0*Bs;
  const T c222 = 3.0*G0;
  const T c223 = 21.0*c115;
  const T c224 = 18.0*c148;
  const T c225 = Bs*c1;
  const T c226 = c115*d;
  const T c227 = c148*d;
  const T c228 = c61*p;
  const T c229 = c225*c5;
  const T c230 = c189*c191;
  const T c231 = pw(Bs, 3)*c1;
  const T c232 = c115*c221;
  const T c233 = c189*c222;
  const T c234 = PA1*d;
  const T c235 = Bs*c151 + Bs*c154 + Bs*c88 + 6.0*c132 + c230 + c231 + c232*c60 + c232 + c233 + c234*c80 + c36*c80;
  const T c236 = c211 + c235;
  const T c237 = 18.0*c115;
  const T c238 = 36.0*c162;
  const T c239 = 36.0*c225;
  const T c240 = 60.0*c1;
  const T c241 = 15.0*c115;
  const T c242 = 24.0*c1;
  const T c243 = Bs*c242;
  const T c244 = 9.0*c115;
  const T c245 = 18.0*c225;
  const T c246 = 18.0*c162;
  const T c247 = Bs*c185;
  K[0] = to_dbl(c2);
  K[1] = to_dbl(c3);
  K[2] = to_dbl(((T(1.0) / T(5040.0)))*c2);
  K[3] = to_dbl(((T(1.0) / T(720.0)))*c2);
  K[4] = to_dbl(((T(1.0) / T(120.0)))*c2);
  K[5] = to_dbl(((T(1.0) / T(24.0)))*c2);
  K[6] = to_dbl(((T(1.0) / T(6.0)))*c2);
  K[7] = to_dbl(d*(((T(1.0) / T(2.0)))*c1 + ((T(3.0) / T(2.0)))*p + 1.0));
  K[8] = to_dbl(c2);
  K[9] = to_dbl(c2);
  K[10] = to_dbl(c3);
  K[11] = to_dbl(c7);
  K[12] = to_dbl(c11);
  K[13] = to_dbl(-(T(1.0) / T(5040.0))*c6);
  K[14] = to_dbl(-(T(1.0) / T(720.0))*c6);
  K[15] = to_dbl(-(T(1.0) / T(120.0))*c6);
  K[16] = to_dbl(-(T(1.0) / T(24.0))*c6);
  K[17] = to_dbl(-(T(1.0) / T(6.0))*c6);
  K[18] = to_dbl(-(T(1.0) / T(2.0))*c6);
  K[19] = to_dbl(c7);
  K[20] = to_dbl(c7);
  K[21] = to_dbl(c11);
  K[22] = to_dbl(c15);
  K[23] = to_dbl(c15);
  K[24] = to_dbl(c19);
  K[25] = to_dbl(c19);
  K[26] = to_dbl(-(T(1.0) / T(5040.0))*c14);
  K[27] = to_dbl(((T(1.0) / T(5040.0)))*c13*(-c20 + c24 - 8.0));
  K[28] = to_dbl(((T(1.0) / T(720.0)))*c13*(c24 - c25 - 7.0));
  K[29] = to_dbl(((T(1.0) / T(120.0)))*c13*(c24 - c26 - 6.0));
  K[30] = to_dbl(((T(1.0) / T(24.0)))*c13*(c24 - c27 - 5.0));
  K[31] = to_dbl(((T(1.0) / T(6.0)))*c13*(c24 - c28 - 4.0));
  K[32] = to_dbl(((T(1.0) / T(2.0)))*c13*(c24 - c29 - 3.0));
  K[33] = to_dbl(c15);
  K[34] = to_dbl(c13*(-c0 + c24 - 2.0));
  K[35] = to_dbl(c30);
  K[36] = to_dbl(c30);
  K[37] = to_dbl(c34);
  K[38] = to_dbl(c34);
  K[39] = to_dbl(((T(1.0) / T(5040.0)))*c34);
  K[40] = to_dbl(((T(1.0) / T(720.0)))*c34);
  K[41] = to_dbl(((T(1.0) / T(120.0)))*c34);
  K[42] = to_dbl(((T(1.0) / T(24.0)))*c34);
  K[43] = to_dbl(((T(1.0) / T(6.0)))*c34);
  K[44] = to_dbl(((T(1.0) / T(2.0)))*c34);
  K[45] = to_dbl(c34);
  K[46] = to_dbl(c34);
  K[47] = to_dbl(c34);
  K[48] = to_dbl(c35);
  K[49] = to_dbl(c35);
  K[50] = to_dbl(c41);
  K[51] = to_dbl(c41);
  K[52] = to_dbl(((T(1.0) / T(5040.0)))*c35);
  K[53] = to_dbl(-c45*(c42 + c44));
  K[54] = to_dbl(-(T(1.0) / T(720.0))*c40*(c44 + c46));
  K[55] = to_dbl(-(T(1.0) / T(120.0))*c40*(c44 + c48));
  K[56] = to_dbl(-(T(1.0) / T(24.0))*c40*(-5.0*c17 + c44));
  K[57] = to_dbl(-(T(1.0) / T(6.0))*c40*(c44 + c50));
  K[58] = to_dbl(-(T(1.0) / T(2.0))*c40*(c44 - c51));
  K[59] = to_dbl(c35);
  K[60] = to_dbl(-c40*(c44 + c53));
  K[61] = to_dbl(c54);
  K[62] = to_dbl(c54);
  K[63] = to_dbl(c40);
  K[64] = to_dbl(c40);
  K[65] = to_dbl(c57);
  K[66] = to_dbl(c57);
  K[67] = to_dbl(c69);
  K[68] = to_dbl(c69);
  K[69] = to_dbl(c45);
  K[70] = to_dbl(-c45*(2.0*c5 + c70 + c71 - 9.0));
  K[71] = to_dbl(((T(1.0) / T(5040.0)))*c68*(-16.0*c17 - 16.0*c21 + c76 - 16.0*c8 + 56.0*p));
  K[72] = to_dbl(((T(1.0) / T(720.0)))*c68*(-14.0*c17 - 14.0*c21 - c5*c77 + c76 + 42.0*p));
  K[73] = to_dbl(((T(1.0) / T(120.0)))*c68*(-12.0*c17 - 12.0*c21 - c5*c79 + c76 + c78));
  K[74] = to_dbl(((T(1.0) / T(24.0)))*c68*(-10.0*c17 - 10.0*c21 + c76 - 10.0*c8 + 20.0*p));
  K[75] = to_dbl(((T(1.0) / T(6.0)))*c68*(-8.0*c21 - c25*c5 + c42 + c76 + c79));
  K[76] = to_dbl(c40);
  K[77] = to_dbl(((T(1.0) / T(2.0)))*c68*(-c27*c5 + c27 + c48 + c76 - c80));
  K[78] = to_dbl(-c56*(Bs + c55));
  K[79] = to_dbl(c68*(c16 - 4.0*c21 - c29*c5 + c50 + c76));
  K[80] = to_dbl(c81);
  K[81] = to_dbl(c81);
  K[82] = to_dbl(c86);
  K[83] = to_dbl(c92);
  K[84] = to_dbl(-(T(1.0) / T(5040.0))*c85);
  K[85] = to_dbl(-(T(1.0) / T(720.0))*c85);
  K[86] = to_dbl(-(T(1.0) / T(120.0))*c85);
  K[87] = to_dbl(-(T(1.0) / T(24.0))*c85);
  K[88] = to_dbl(-(T(1.0) / T(6.0))*c85);
  K[89] = to_dbl(-(T(1.0) / T(2.0))*c85);
  K[90] = to_dbl(c86);
  K[91] = to_dbl(c86);
  K[92] = to_dbl(c92);
  K[93] = to_dbl(c94);
  K[94] = to_dbl(c94);
  K[95] = to_dbl(c100*c102);
  K[96] = to_dbl(c101*c114);
  K[97] = to_dbl(-(T(1.0) / T(5040.0))*c93);
  K[98] = to_dbl(c121*(-7.0*c115 + c117 + c120));
  K[99] = to_dbl(((T(1.0) / T(720.0)))*c102*(-G1*c28 + c120 + c123));
  K[100] = to_dbl(((T(1.0) / T(120.0)))*c102*(-5.0*c115 + c120 + c124));
  K[101] = to_dbl(((T(1.0) / T(24.0)))*c102*(c119 + c126 + c84 + c91 + c98));
  K[102] = to_dbl(((T(1.0) / T(6.0)))*c102*(c110 + c120 + c128));
  K[103] = to_dbl(((T(1.0) / T(2.0)))*c102*(c120 + c129 + c131));
  K[104] = to_dbl(c94);
  K[105] = to_dbl(c102*(-c115 + c120));
  K[106] = to_dbl(c102*(c100 + c119));
  K[107] = to_dbl(c101*(c114 + c132 + c21*c36));
  K[108] = to_dbl(c134);
  K[109] = to_dbl(c134);
  K[110] = to_dbl(c136);
  K[111] = to_dbl(c136);
  K[112] = to_dbl(-c144*c146);
  K[113] = to_dbl(-c145*c159);
  K[114] = to_dbl(-(T(1.0) / T(5040.0))*c133);
  K[115] = to_dbl(c121*(2.0*G1 + c160 + c58 + c62));
  K[116] = to_dbl(-(T(1.0) / T(5040.0))*c146*(-c161*d - c161 + c163 - 14.0*c164 + c165 + c170 + c172));
  K[117] = to_dbl(-(T(1.0) / T(720.0))*c146*(-10.0*c148 + c170 + c174 + c175 - 10.0*c176 + c178 - 10.0*c95));
  K[118] = to_dbl(-(T(1.0) / T(120.0))*c146*(-G1*c25 - PA1*c25 + 20.0*c162 - 10.0*c164 + c170 - c179*d - c179 - c25*c36));
  K[119] = to_dbl(-(T(1.0) / T(24.0))*c146*(-8.0*c164 + c170 + c177 - c180*d - c180 + c182 + c184));
  K[120] = to_dbl(-(T(1.0) / T(6.0))*c146*(c124 - c139 + c150 + c170 + c187 + c188));
  K[121] = to_dbl(c134);
  K[122] = to_dbl(-(T(1.0) / T(2.0))*c146*(-c106 + c107 + c113 - c125*d + c126 - c137 + c142 - 4.0*c164 + c169));
  K[123] = to_dbl(c135*c44);
  K[124] = to_dbl(-c146*(-c130*d + c131 - c162*c70 + c170));
  K[125] = to_dbl(-c146*(c144 + c169));
  K[126] = to_dbl(-c145*(c159 + c166*d + c167*c60 + c168 + c189*c5 + c61*c74));
  K[127] = to_dbl(c190);
  K[128] = to_dbl(c190);
  K[129] = to_dbl(c192);
  K[130] = to_dbl(c192);
  K[131] = to_dbl(c194);
  K[132] = to_dbl(c194);
  K[133] = to_dbl(c212*c214);
  K[134] = to_dbl(c213*c220);
  K[135] = to_dbl(-c121);
  K[136] = to_dbl(c121*(c191 + c221 + c222 - 7.0));
  K[137] = to_dbl(-(T(1.0) / T(1680.0))*c146*(c160 - 7.0*c21 - c26*c5 + c75 + c77));
  K[138] = to_dbl(((T(1.0) / T(5040.0)))*c214*(-Bs*c163 - 210.0*c1 + 126.0*c149 + 126.0*c162 - 36.0*c176 - 21.0*c189 - c223*c60 - c223 - c224 + 126.0*c225 - 42.0*c226 - 36.0*c227 - 18.0*c228 - 42.0*c229 + c236 - 18.0*c95));
  K[139] = to_dbl(((T(1.0) / T(720.0)))*c214*(-Bs*c238 - 120.0*c1 - 15.0*c148 + 90.0*c149 + 90.0*c162 - 18.0*c189 + 90.0*c225 - 36.0*c226 - 15.0*c228 - c234*c78 + c236 - c237*c60 - c237 - c239*c5 - c36*c78 - 15.0*c95));
  K[140] = to_dbl(((T(1.0) / T(120.0)))*c214*(-Bs*c175 + Bs*c240 + G0*c240 + c172 - 24.0*c176 - 15.0*c189 - 30.0*c226 - 24.0*c227 - 30.0*c229 + c236 + c240*c5 - c240 - c241*c60 - c241 - c61*c79));
  K[141] = to_dbl(((T(1.0) / T(24.0)))*c214*(-G0*c243 + 36.0*c149 - c173*c60 + c174 - 18.0*c176 - 12.0*c189 - c20*c61 - c203 - c205 - c224*d - 24.0*c226 + c236 + c238 + c239 - c242 - c243*c5));
  K[142] = to_dbl(c190);
  K[143] = to_dbl(((T(1.0) / T(6.0)))*c214*(-Bs*c246 + 18.0*c149 + c165 - c171*d + c184 - c185 - 9.0*c189 + c196 + c197 + c198 + c199 + c200 + c201 + c202 + c204 + c206 + c210 + c219 + c235 - c237*d - c244*c60 - c244 - c245*c5 + c245 + c246));
  K[144] = to_dbl(3.0*c102*c23);
  K[145] = to_dbl(((T(1.0) / T(2.0)))*c214*(-c122*c60 + c140 - c151 + c158 + c178 + c182 - c185*c73 + c188 - c195 + c209 + c216 - 12.0*c229 + c235 + c247 + c84 + c90));
  K[146] = to_dbl(-c193*c76);
  K[147] = to_dbl(c214*(-c127*c60 + c128 + c187 - 3.0*c189 + c236 - c247*c5));
  K[148] = to_dbl(c214*(c212 + c235));
  K[149] = to_dbl(c213*(Bs*c104 + Bs*c217 + Bs*c218 + 6.0*c118*c60 + c153*c80 + c155*c232 + c220 + c221*c226 + c230 + c231*d + c233*c60 + c61*c80));
}

template <class T>
inline void gengamma2_coef_deriv4(const T& d, const T& p, const T* PAv, const T* Gv, const T& Bs, double* K) {
  (void) d; (void) p; (void) PAv; (void) Gv; (void) Bs;
  const T PA1 = PAv[1];
  const T PA2 = PAv[2];
  const T PA3 = PAv[3];
  const T G0 = Gv[0];
  const T G1 = Gv[1];
  const T G2 = Gv[2];
  const T G3 = Gv[3];
  const T c0 = 11.0*p;
  const T c1 = pw(p, 3);
  const T c2 = pw(p, 2);
  const T c3 = 6.0*c2;
  const T c4 = d*(c0 + c1 + c3 + 6.0);
  const T c5 = -c4;
  const T c6 = 6.0*p;
  const T c7 = d*p;
  const T c8 = -c7*(c2 + c6 + 11.0);
  const T c9 = 3.0*c2;
  const T c10 = 3.0*p;
  const T c11 = c10 + c2 + 2.0;
  const T c12 = G0*d;
  const T c13 = c11*c12;
  const T c14 = 2.0*G0;
  const T c15 = c14*d;
  const T c16 = G0*c2;
  const T c17 = c16*d;
  const T c18 = c10*c12 + c15 + c17;
  const T c19 = c18 - 2.0;
  const T c20 = pw(p, -1);
  const T c21 = c20*d;
  const T c22 = c11*c21;
  const T c23 = G0*c10;
  const T c24 = c16 + c23;
  const T c25 = c19 + c24 - c6 - c9;
  const T c26 = -c21*(c14 + c25);
  const T c27 = 10.0*c2;
  const T c28 = Bs*c2;
  const T c29 = Bs*c10;
  const T c30 = 2.0*Bs;
  const T c31 = c14 + c30;
  const T c32 = c28 + c29 + c31;
  const T c33 = c18 + c24 + c32;
  const T c34 = 24.0*p;
  const T c35 = 9.0*c2;
  const T c36 = 21.0*p;
  const T c37 = 8.0*c2;
  const T c38 = 18.0*p;
  const T c39 = 7.0*c2;
  const T c40 = 15.0*p;
  const T c41 = -c3;
  const T c42 = 12.0*p;
  const T c43 = 5.0*c2;
  const T c44 = 9.0*p;
  const T c45 = 4.0*c2;
  const T c46 = -c21*(c25 + c32);
  const T c47 = G1*p;
  const T c48 = pw(G0, 2);
  const T c49 = c2*c48;
  const T c50 = c48*p;
  const T c51 = G1 + c50;
  const T c52 = c47 + c49 + c51;
  const T c53 = c21*c52;
  const T c54 = -c53;
  const T c55 = p + 1.0;
  const T c56 = c12*c20*c55;
  const T c57 = -c56;
  const T c58 = G0*p;
  const T c59 = -c58;
  const T c60 = c14*c2;
  const T c61 = c47*d;
  const T c62 = c49*d;
  const T c63 = G1*d;
  const T c64 = c50*d;
  const T c65 = c63 + c64;
  const T c66 = PA1*p;
  const T c67 = PA1 + c52 + c66;
  const T c68 = c61 + c62 + c65 + c67;
  const T c69 = c59 - c60 + c68;
  const T c70 = pw(c2, -1);
  const T c71 = c70*d;
  const T c72 = c69*c71;
  const T c73 = 9.0*c16;
  const T c74 = 8.0*c58;
  const T c75 = -c74;
  const T c76 = Bs*c58;
  const T c77 = Bs*c16;
  const T c78 = c76 + c77;
  const T c79 = c68 + c78;
  const T c80 = ((T(1.0) / T(5040.0)))*c71;
  const T c81 = -7.0*c58;
  const T c82 = 8.0*c16;
  const T c83 = -c82;
  const T c84 = G0*c6;
  const T c85 = -c84;
  const T c86 = G0*c3;
  const T c87 = -c86;
  const T c88 = 4.0*c58;
  const T c89 = -c88;
  const T c90 = 4.0*c16;
  const T c91 = -c90;
  const T c92 = c14*p;
  const T c93 = -c92;
  const T c94 = G0*c9;
  const T c95 = c71*(c69 + c78);
  const T c96 = c55*c71;
  const T c97 = -c96;
  const T c98 = G0 + c12;
  const T c99 = c98 - 1.0;
  const T c100 = 2.0*p;
  const T c101 = c58*d;
  const T c102 = -c100 + c101 + c58;
  const T c103 = 2.0*c71;
  const T c104 = c103*(c102 + c99);
  const T c105 = 2.0*c63 + 2.0*c64;
  const T c106 = c14*c7;
  const T c107 = -c106 + c93;
  const T c108 = 2.0*c2;
  const T c109 = 4.0*c17;
  const T c110 = 2.0*c61;
  const T c111 = pw(d, 2);
  const T c112 = c111*c47;
  const T c113 = c111*c49;
  const T c114 = 2.0*c66;
  const T c115 = 2.0*c62;
  const T c116 = G1*c111;
  const T c117 = 2.0*PA1;
  const T c118 = c111*c50 + c116 + c117*d;
  const T c119 = c110 + c112 + c113 + c114*d + c115 + c118 + c67;
  const T c120 = c108 - c109 + c119 + c91;
  const T c121 = d/c1;
  const T c122 = -c121*(c105 + c107 + c120);
  const T c123 = c30*p;
  const T c124 = 72.0*c2;
  const T c125 = 18.0*c28;
  const T c126 = 18.0*c16;
  const T c127 = Bs*p;
  const T c128 = 16.0*c58;
  const T c129 = 18.0*c17;
  const T c130 = 16.0*c101;
  const T c131 = pw(Bs, 2);
  const T c132 = c131*p;
  const T c133 = c127*c15;
  const T c134 = c131*c2;
  const T c135 = c14*c28;
  const T c136 = c105 + c127*c14;
  const T c137 = c132 + c133 + c134 + c135 + c136 + c15*c28;
  const T c138 = c119 + c137;
  const T c139 = ((T(1.0) / T(5040.0)))*c121;
  const T c140 = 16.0*c16;
  const T c141 = -14.0*c101 - 14.0*c127 - 14.0*c58;
  const T c142 = ((T(1.0) / T(720.0)))*c121;
  const T c143 = 30.0*p;
  const T c144 = 14.0*c16;
  const T c145 = Bs*c42;
  const T c146 = G0*c42;
  const T c147 = ((T(1.0) / T(120.0)))*c121;
  const T c148 = 20.0*p;
  const T c149 = 30.0*c2;
  const T c150 = 12.0*c28;
  const T c151 = 12.0*c16;
  const T c152 = 12.0*c17;
  const T c153 = ((T(1.0) / T(24.0)))*c121;
  const T c154 = 20.0*c2;
  const T c155 = Bs*c27;
  const T c156 = ((T(1.0) / T(6.0)))*c121;
  const T c157 = 12.0*c2;
  const T c158 = 8.0*c28;
  const T c159 = Bs*c6;
  const T c160 = ((T(1.0) / T(2.0)))*c121;
  const T c161 = Bs + c99;
  const T c162 = Bs*c3;
  const T c163 = 4.0*c127;
  const T c164 = c12*c3;
  const T c165 = c88*d;
  const T c166 = 4.0*c28;
  const T c167 = c107 - c123;
  const T c168 = -c121*(c120 + c137 - c166 + c167);
  const T c169 = G1*c23;
  const T c170 = pw(G0, 3);
  const T c171 = c170*c2;
  const T c172 = G2 + c169 + c171;
  const T c173 = c172*c71;
  const T c174 = c51*c71;
  const T c175 = G2*d;
  const T c176 = PA2 + c172;
  const T c177 = c14*c66 + c171*d + c175 + c176 + c23*c63;
  const T c178 = c177 - c49;
  const T c179 = -c121*c178;
  const T c180 = c37*c48;
  const T c181 = -c180;
  const T c182 = Bs*c47 + c28*c48;
  const T c183 = c177 + c182;
  const T c184 = c39*c48;
  const T c185 = -c184;
  const T c186 = G1*c6;
  const T c187 = -c186;
  const T c188 = c3*c48;
  const T c189 = -c188;
  const T c190 = 5.0*c47;
  const T c191 = -c190;
  const T c192 = c43*c48;
  const T c193 = -4.0*c47;
  const T c194 = c45*c48;
  const T c195 = -c194;
  const T c196 = G1*c10;
  const T c197 = -c196;
  const T c198 = c48*c9;
  const T c199 = -c198;
  const T c200 = 2.0*c47;
  const T c201 = -c200;
  const T c202 = -2.0*c49;
  const T c203 = -c121*(c178 + c182);
  const T c204 = c12*c70;
  const T c205 = PA1 + c51;
  const T c206 = c205 + c65;
  const T c207 = c206 + c59;
  const T c208 = 2.0*c121;
  const T c209 = -c207*c208;
  const T c210 = c108*c170;
  const T c211 = 2.0*c175 + c210*d + c63*c84;
  const T c212 = G2*c111;
  const T c213 = 2.0*PA2;
  const T c214 = 4.0*c66;
  const T c215 = PA1*c23 + c176;
  const T c216 = c111*c171 + c116*c23 + c12*c214 + c212 + c213*d + c215;
  const T c217 = c211 + c216;
  const T c218 = -c115 + c202;
  const T c219 = pw(p, 4);
  const T c220 = pw(c219, -1);
  const T c221 = c220*d;
  const T c222 = c221*(c217 + c218);
  const T c223 = G0*c44;
  const T c224 = c136 - c223;
  const T c225 = c211 + c30*c47 + c30*c49;
  const T c226 = c131*c16 + c216 + c30*c61 + c30*c62 + c30*c66;
  const T c227 = c225 + c226;
  const T c228 = -c42*c63;
  const T c229 = G1*c42;
  const T c230 = PA1*c42;
  const T c231 = -c229 - c230;
  const T c232 = 12.0*c49;
  const T c233 = -c232;
  const T c234 = 30.0*c16;
  const T c235 = 12.0*c77;
  const T c236 = -c235 - 12.0*c62;
  const T c237 = c27*c48;
  const T c238 = 20.0*c16;
  const T c239 = c6*c63;
  const T c240 = -c239;
  const T c241 = PA1*c6;
  const T c242 = c187 - c241;
  const T c243 = c188*d;
  const T c244 = -c243;
  const T c245 = Bs*c86;
  const T c246 = -c245;
  const T c247 = c189 + c86;
  const T c248 = 4.0*Bs;
  const T c249 = c16*c248;
  const T c250 = c221*(-c135 + c218 + c227);
  const T c251 = 3.0*c121;
  const T c252 = -c251*c99;
  const T c253 = c118 + c205;
  const T c254 = c105 + c253;
  const T c255 = 3.0*c221;
  const T c256 = c255*(c107 + c254);
  const T c257 = c111*c198;
  const T c258 = pw(d, 3);
  const T c259 = G2*c258;
  const T c260 = 3.0*PA2;
  const T c261 = c170*c9;
  const T c262 = PA1*c84;
  const T c263 = PA1*c44;
  const T c264 = c111*c260 + c111*c261 + c111*c262 + c116*c223 + c12*c263 + c169*c258 + c171*c258 + 3.0*c175 + 3.0*c212 + c215 + c223*c63 + c259 + c260*d + c261*d;
  const T c265 = c199 + c244 - c257 + c264;
  const T c266 = pw(p, -5);
  const T c267 = c266*d;
  const T c268 = -c265*c267;
  const T c269 = 3.0*Bs;
  const T c270 = 3.0*G0;
  const T c271 = c132 + c133 + c253;
  const T c272 = 24.0*c49;
  const T c273 = c66*d;
  const T c274 = 24.0*c113;
  const T c275 = Bs*c17;
  const T c276 = pw(Bs, 3);
  const T c277 = c131*c9;
  const T c278 = PA1*d;
  const T c279 = c48*d;
  const T c280 = Bs*c198 + Bs*c257 + G1*c29 + PA1*c29 + c116*c29 + c12*c277 + c131*c94 + c159*c278 + c159*c63 + c162*c279 + c2*c276;
  const T c281 = c264 + c280;
  const T c282 = ((T(1.0) / T(5040.0)))*c267;
  const T c283 = PA1*c38;
  const T c284 = c116*c38;
  const T c285 = -21.0*c113 - 21.0*c134 - 42.0*c275 - 21.0*c49 - 42.0*c62 - 42.0*c77;
  const T c286 = ((T(1.0) / T(720.0)))*c267;
  const T c287 = 120.0*c2;
  const T c288 = 18.0*c49;
  const T c289 = 90.0*c16;
  const T c290 = 36.0*c77;
  const T c291 = 36.0*c17;
  const T c292 = Bs*c291;
  const T c293 = PA1*c40;
  const T c294 = -G1*c40 - c116*c40 - 30.0*c273 - c293 - 30.0*c61;
  const T c295 = ((T(1.0) / T(120.0)))*c267;
  const T c296 = 60.0*c2;
  const T c297 = 60.0*c28;
  const T c298 = 60.0*c16;
  const T c299 = c34*c63;
  const T c300 = c34*d;
  const T c301 = 15.0*c113;
  const T c302 = 60.0*c17;
  const T c303 = ((T(1.0) / T(24.0)))*c267;
  const T c304 = 24.0*c2;
  const T c305 = 12.0*c134;
  const T c306 = G1*c44;
  const T c307 = 36.0*c28;
  const T c308 = 36.0*c16;
  const T c309 = 24.0*c77;
  const T c310 = 24.0*c62;
  const T c311 = c38*c63;
  const T c312 = ((T(1.0) / T(6.0)))*c267;
  const T c313 = c35*c48;
  const T c314 = 18.0*c77;
  const T c315 = c111*c313;
  const T c316 = ((T(1.0) / T(2.0)))*c267;
  const T c317 = c136 + c271;
  const T c318 = c131*c3;
  const T c319 = PA1*c10;
  const T c320 = c241*d;
  const T c321 = c111*c188;
  const T c322 = c10*c116;
  const T c323 = Bs*c152;
  const T c324 = -c267*(-c12*c162 + c246 + c265 - c277 + c280);
  const T c325 = G2*c88;
  const T c326 = pw(G1, 2);
  const T c327 = c10*c326;
  const T c328 = pw(G0, 4)*c1;
  const T c329 = G1*c188;
  const T c330 = G3 + c325 + c327 + c328 + c329;
  const T c331 = c221*c330;
  const T c332 = -c331;
  const T c333 = G2*p;
  const T c334 = 4.0*c333;
  const T c335 = PA3 - c334;
  const T c336 = c328*d;
  const T c337 = c175*c88 + c188*c63 + c327*d + c336;
  const T c338 = c335 + c337;
  const T c339 = -c220*c338;
  const T c340 = c172*c221;
  const T c341 = -c340;
  const T c342 = G1*c94;
  const T c343 = G2*c100;
  const T c344 = PA3 + c330;
  const T c345 = c343 + c344;
  const T c346 = G3*d;
  const T c347 = PA1*c198;
  const T c348 = PA1*c196 + PA2*c23 + c346 + c347;
  const T c349 = c337 + c348;
  const T c350 = c342 + c345 + c349;
  const T c351 = PA3*d;
  const T c352 = G2*c10;
  const T c353 = -c352;
  const T c354 = c63*c94;
  const T c355 = PA2*c10;
  const T c356 = G1*c3;
  const T c357 = c337 - c356;
  const T c358 = 4.0*p;
  const T c359 = c175*c358;
  const T c360 = -c359;
  const T c361 = c111*c328;
  const T c362 = c111*c327 + c116*c188 + c212*c88 + c361;
  const T c363 = c360 + c362;
  const T c364 = c12*c355 + c319*c63 + c347*d + c351 + c353 + c354 + c357 + c363;
  const T c365 = c1*c170;
  const T c366 = -G1*c126;
  const T c367 = Bs*c365;
  const T c368 = Bs*c342 + G2*c127 + c367;
  const T c369 = c344 + c349 + c368;
  const T c370 = 6.0*c365;
  const T c371 = -c370;
  const T c372 = G1*c16;
  const T c373 = G1*c151;
  const T c374 = -c373;
  const T c375 = 4.0*c365;
  const T c376 = -c375;
  const T c377 = -c343;
  const T c378 = G1*c73;
  const T c379 = -c378;
  const T c380 = -c333;
  const T c381 = 3.0*c365;
  const T c382 = -c381;
  const T c383 = G1*c86;
  const T c384 = -c383;
  const T c385 = 2.0*c365;
  const T c386 = -c385;
  const T c387 = c367*d;
  const T c388 = c221*c51;
  const T c389 = -c388;
  const T c390 = c177 + c47;
  const T c391 = 2.0*c267;
  const T c392 = c390*c391;
  const T c393 = PA1*c90;
  const T c394 = c63*c86;
  const T c395 = PA2*p;
  const T c396 = 4.0*c395;
  const T c397 = PA2*c88;
  const T c398 = c344 + c397;
  const T c399 = c334 + c396 + c398;
  const T c400 = G3*c111;
  const T c401 = pw(PA1, 2);
  const T c402 = c100*c401;
  const T c403 = PA1*c192;
  const T c404 = c326*c6;
  const T c405 = PA2*c6;
  const T c406 = PA1*c243 + c241*c63;
  const T c407 = PA1*c190 + c12*c405 + c175*c74 + c232*c63 + 2.0*c336 + 2.0*c346 + 2.0*c351 + c400 + c402 + c403 + c404*d + c406;
  const T c408 = c362 + c407;
  const T c409 = G1*c108 + c408;
  const T c410 = c359 + c383 + c393 + c394 + c399 + c409;
  const T c411 = pw(p, -6);
  const T c412 = c411*d;
  const T c413 = PA3*c111;
  const T c414 = 4.0*G0;
  const T c415 = c1*c414;
  const T c416 = 5.0*c61;
  const T c417 = PA2*c111;
  const T c418 = PA2*c100;
  const T c419 = c377 - c418;
  const T c420 = c175*c6;
  const T c421 = -c420;
  const T c422 = PA2*c165;
  const T c423 = c421 + c422;
  const T c424 = c258*c328;
  const T c425 = c258*c327 + c258*c329 + c259*c88 + c424;
  const T c426 = -G1*c37 - PA1*c108 + PA1*c109 + PA1*c321 + PA1*c416 + c111*c404 + c116*c232 + c116*c241 + c116*c86 - c157*c63 - c212*c358 + c212*c74 + c337 + 2.0*c361 + c394 + c402*d + c403*d + c413 - c415 + c417*c84 + c419 + c423 + c425;
  const T c427 = G1*c84;
  const T c428 = 2.0*G2 + c210 + c213 + c427;
  const T c429 = c185 + c191 + c225;
  const T c430 = 14.0*c365;
  const T c431 = c1*c48;
  const T c432 = 42.0*c431;
  const T c433 = G1*c308;
  const T c434 = -c433;
  const T c435 = PA1*c16;
  const T c436 = 24.0*c435;
  const T c437 = -c436;
  const T c438 = Bs*c431;
  const T c439 = c175*p;
  const T c440 = c131*c431;
  const T c441 = c30*c395;
  const T c442 = c30*c365;
  const T c443 = c123*c175 + c245*c63 + c442*d;
  const T c444 = G1*c134 + G1*c245 + PA1*c249 + c30*c333 + c440 + c441 + c442 + c443;
  const T c445 = c408 + c444;
  const T c446 = c398 + c445;
  const T c447 = G1*c150;
  const T c448 = c308*c63;
  const T c449 = -c447 - c448;
  const T c450 = 12.0*c365;
  const T c451 = -c450;
  const T c452 = 30.0*c431;
  const T c453 = 12.0*c431;
  const T c454 = Bs*c453;
  const T c455 = -c450*d - c454;
  const T c456 = 8.0*c333;
  const T c457 = 8.0*c395;
  const T c458 = -c456 - c457;
  const T c459 = 10.0*c365;
  const T c460 = c16*c63;
  const T c461 = G2*c6;
  const T c462 = -c405 - c461;
  const T c463 = G1*c157;
  const T c464 = 24.0*c372;
  const T c465 = c463 - c464;
  const T c466 = 8.0*c365;
  const T c467 = G1*c162;
  const T c468 = c126*c63;
  const T c469 = -c467 - c468;
  const T c470 = c356 + c366;
  const T c471 = PA1*c151;
  const T c472 = c453 - c466 - c471;
  const T c473 = c330 + c335 - c396 + c397;
  const T c474 = c100*c175;
  const T c475 = c151*c63;
  const T c476 = c374 + c398;
  const T c477 = 6.0*c431;
  const T c478 = Bs*c477;
  const T c479 = -c370*d - c478;
  const T c480 = c371 + c477;
  const T c481 = c2*c30;
  const T c482 = G1*c481;
  const T c483 = c440*d;
  const T c484 = PA1*c17;
  const T c485 = c12*c220;
  const T c486 = -c485;
  const T c487 = 3.0*c267;
  const T c488 = c206*c487;
  const T c489 = c110 + c114 + c200;
  const T c490 = 3.0*c412;
  const T c491 = -c490*(c217 + c489);
  const T c492 = c212*c6;
  const T c493 = PA1*c73;
  const T c494 = c116*c73;
  const T c495 = PA1*c152;
  const T c496 = c175*c42;
  const T c497 = PA2*c42;
  const T c498 = c497*d;
  const T c499 = c496 + c498;
  const T c500 = PA1*c3;
  const T c501 = c326*c44;
  const T c502 = c111*c501 + c116*c288 + c146*c212 + 3.0*c361;
  const T c503 = G3*c258;
  const T c504 = c401*c6;
  const T c505 = G2*c42;
  const T c506 = c10*c401;
  const T c507 = PA1*c186 + PA1*c188 + c506;
  const T c508 = c398 + c507;
  const T c509 = PA1*c315 + 15.0*PA1*c62 + c116*c263 + c12*c497 + c12*c505 + c223*c417 + c288*c63 + c293*c63 + 3.0*c336 + 3.0*c346 + 3.0*c351 + 3.0*c400 + 3.0*c413 + c425 + c501*d + c503 + c504*d + c508;
  const T c510 = c502 + c509;
  const T c511 = c3*c63 + c500 + c510;
  const T c512 = c356 + c378 + c405 + c461 + c468 + c492 + c493 + c494 + c495 + c499 + c511;
  const T c513 = pw(p, -7);
  const T c514 = c513*d;
  const T c515 = PA3*c258;
  const T c516 = 18.0*c2;
  const T c517 = c1*c12;
  const T c518 = c63*c73;
  const T c519 = c258*c378;
  const T c520 = PA2*c258;
  const T c521 = PA1*c258;
  const T c522 = c116*c126;
  const T c523 = G1*c258;
  const T c524 = c304*c63;
  const T c525 = c406 + c506*d - c524;
  const T c526 = -c212*c44 + c502;
  const T c527 = c259*c358;
  const T c528 = pw(d, 4);
  const T c529 = c328*c528;
  const T c530 = c325*c528 + c327*c528 + c329*c528 + c529;
  const T c531 = -c527 + c530;
  const T c532 = G0*c1;
  const T c533 = 6.0*c532;
  const T c534 = c405*d;
  const T c535 = -c533 - c534;
  const T c536 = PA1*c301 + c111*c471 + c111*c504 + c116*c293 - c116*c516 + c146*c259 + c146*c417 + c223*c520 + c258*c501 + c288*c523 + c306*c521 + c313*c521 + c357 + c380 - c395 + c423 + 3.0*c424 + 9.0*c484 - c500*d - c500 + c515 - 12.0*c517 + c518 + c519 + c522 + c525 + c526 + c531 + c535;
  const T c537 = c10*c48;
  const T c538 = 3.0*G1 + 3.0*PA1 + c537;
  const T c539 = 21.0*c365;
  const T c540 = 90.0*c2;
  const T c541 = 126.0*c431;
  const T c542 = Bs*c432;
  const T c543 = c365*d;
  const T c544 = PA2*d;
  const T c545 = c131*c532;
  const T c546 = Bs*c532;
  const T c547 = -72.0*c484;
  const T c548 = c116*c16;
  const T c549 = c269*c365;
  const T c550 = PA1*c277;
  const T c551 = 3.0*c440;
  const T c552 = 6.0*d;
  const T c553 = PA1*c323;
  const T c554 = c277*c63 + 3.0*c483;
  const T c555 = Bs*c378 + Bs*c493 + Bs*c494 + G1*c277 + G2*c29 + PA2*c29 + c111*c549 + c159*c175 + c159*c544 + c212*c29 + c276*c532 + c314*c63 + c367*c552 + c549 + c550 + c551 + c553 + c554;
  const T c556 = c510 + c555;
  const T c557 = -G1*c307 - PA1*c307 - 108.0*c460;
  const T c558 = c212*c42;
  const T c559 = -c558;
  const T c560 = 120.0*c1;
  const T c561 = 18.0*c365;
  const T c562 = 90.0*c431;
  const T c563 = 30.0*c28;
  const T c564 = 36.0*c431;
  const T c565 = Bs*c564;
  const T c566 = -36.0*c543 - c565;
  const T c567 = -c497 - c505;
  const T c568 = c175*c34;
  const T c569 = PA2*c300;
  const T c570 = -c568 - c569;
  const T c571 = 60.0*c532;
  const T c572 = 15.0*c365;
  const T c573 = 60.0*c431;
  const T c574 = Bs*c452;
  const T c575 = 36.0*c2;
  const T c576 = c116*c308;
  const T c577 = 24.0*c28;
  const T c578 = c577*c63;
  const T c579 = -c576 - c578;
  const T c580 = PA1*c577;
  const T c581 = -G1*c577 - 72.0*c460 - c580;
  const T c582 = PA1*c308;
  const T c583 = G1*c575 + PA1*c575 + c434 - c582;
  const T c584 = 24.0*c1;
  const T c585 = Bs*c584;
  const T c586 = c48*c585;
  const T c587 = -c111*c450 - 12.0*c545 - c586*d;
  const T c588 = 24.0*c365;
  const T c589 = PA1*c291;
  const T c590 = -c586 - c588*d - c589;
  const T c591 = G0*c584;
  const T c592 = c451 - c496 - c498 + c564 - c591;
  const T c593 = 9.0*c365;
  const T c594 = 18.0*c431;
  const T c595 = Bs*c594;
  const T c596 = c10*c212;
  const T c597 = c150*c63;
  const T c598 = PA1*c150;
  const T c599 = c449 - c598;
  const T c600 = PA1*c162;
  const T c601 = c162*c63;
  const T c602 = c1*c131;
  const T c603 = c1*c276;
  const T c604 = c111*c367;
  const T c605 = PA1*c111;
  const T c606 = -c267;
  const T c607 = 4.0*c12;
  const T c608 = c266*c607*(d + 1.0);
  const T c609 = 6.0*c412;
  const T c610 = -c254*c609;
  const T c611 = c196 + c239 + c264 + c319 + c320 + c322;
  const T c612 = 4.0*c514;
  const T c613 = c611*c612;
  const T c614 = c259*p;
  const T c615 = 8.0*c614;
  const T c616 = c212*c34;
  const T c617 = c34*c417;
  const T c618 = c258*c373;
  const T c619 = c111*c436;
  const T c620 = c401*c42;
  const T c621 = c111*c620;
  const T c622 = c576 + c621;
  const T c623 = G3*c528;
  const T c624 = 4.0*c346;
  const T c625 = 4.0*c503;
  const T c626 = 4.0*c351;
  const T c627 = 4.0*c515;
  const T c628 = 6.0*c400;
  const T c629 = 6.0*c413;
  const T c630 = 4.0*c336;
  const T c631 = 4.0*c424;
  const T c632 = 6.0*c361;
  const T c633 = c42*d;
  const T c634 = c326*c633;
  const T c635 = c401*c633;
  const T c636 = c326*c42;
  const T c637 = c258*c636;
  const T c638 = c326*c38;
  const T c639 = c111*c638;
  const T c640 = c146*c520;
  const T c641 = c229*c521;
  const T c642 = c232*c521;
  const T c643 = c128*c175;
  const T c644 = c128*c259;
  const T c645 = PA2*c130;
  const T c646 = G0*c616;
  const T c647 = G0*c617;
  const T c648 = PA1*c299;
  const T c649 = c272*c63;
  const T c650 = c272*c523;
  const T c651 = PA1*c310;
  const T c652 = 30.0*PA1;
  const T c653 = c112*c652;
  const T c654 = c113*c652;
  const T c655 = 36.0*c49;
  const T c656 = c116*c655;
  const T c657 = c530 + c623 + c624 + c625 + c626 + c627 + c628 + c629 + c630 + c631 + c632 + c634 + c635 + c637 + c639 + c640 + c641 + c642 + c643 + c644 + c645 + c646 + c647 + c648 + c649 + c650 + c651 + c653 + c654 + c656;
  const T c658 = c508 + c657;
  const T c659 = PA1*c157;
  const T c660 = c278*c304;
  const T c661 = c116*c157 + c524 + c659 + c660;
  const T c662 = c373 + c448 + c456 + c457 + c463 + c471 + c568 + c569 + c589 + c615 + c616 + c617 + c618 + c619 + c622 + c658 + c661;
  const T c663 = pw(p, -8);
  const T c664 = c663*d;
  const T c665 = pw(d, 5);
  const T c666 = 48.0*c116;
  const T c667 = -c259*c42;
  const T c668 = PA1*c528;
  const T c669 = G0*c34;
  const T c670 = G1*c528;
  const T c671 = 30.0*c521;
  const T c672 = c258*c433;
  const T c673 = c12*c584;
  const T c674 = c111*c497;
  const T c675 = c559 - c673 - c674;
  const T c676 = G2*c128*c528 + PA1*c116*c34 + PA1*c274 + PA2*c146*c528 + PA3*c528 + c111*c582 - c111*c591 + c111*c636 - c111*c659 + c116*c272 + c128*c212 + c128*c417 - c2*c666 + 6.0*c219 + c229*c668 + c232*c668 + c258*c436 + c258*c620 + c258*c638 + c259*c669 + c272*c670 - c304*c523 + c325*c665 + c327*c665 + c328*c665 + c329*c665 - c334*c528 + c337 + c360 + 4.0*c361 + c373*c528 - c396*d + c422 + 6.0*c424 + c47*c671 + c475 + c49*c671 + c495 + c520*c669 + c523*c655 + c525 + c528*c636 + 4.0*c529 + c622 - c660 + c667 + c672 + c675;
  const T c677 = c162*c48;
  const T c678 = c170*c3;
  const T c679 = 840.0*c1;
  const T c680 = 28.0*c365;
  const T c681 = 180.0*c2;
  const T c682 = 252.0*c431;
  const T c683 = 84.0*c438;
  const T c684 = c111*c365;
  const T c685 = 72.0*c28;
  const T c686 = 72.0*c372;
  const T c687 = 72.0*c435;
  const T c688 = c212*p;
  const T c689 = 60.0*c395;
  const T c690 = 360.0*c2;
  const T c691 = c431*d;
  const T c692 = c438*d;
  const T c693 = 144.0*c28;
  const T c694 = c111*c435;
  const T c695 = c12*c602;
  const T c696 = Bs*c1;
  const T c697 = c12*c696;
  const T c698 = pw(Bs, 4)*c1;
  const T c699 = G2*c163;
  const T c700 = c248*c365;
  const T c701 = c276*c415;
  const T c702 = 6.0*c440;
  const T c703 = c116*c290 + c145*c212 + c145*c417 + c603*c607 + 12.0*c604;
  const T c704 = G1*c235 + G1*c318 + PA1*c235 + PA1*c292 + PA1*c318 + PA2*c163 + c111*c702 + c116*c318 + c145*c175 + c145*c544 + c163*c259 + c235*c523 + c258*c700 + c278*c305 + c290*c63 + c305*c63 + c309*c605 + 12.0*c387 + 12.0*c483 + c698 + c699 + c700 + c701 + c702 + c703;
  const T c705 = c621 + c704;
  const T c706 = c658 + c705;
  const T c707 = 360.0*c1;
  const T c708 = 180.0*c431;
  const T c709 = 72.0*c431;
  const T c710 = Bs*c709;
  const T c711 = 48.0*c395;
  const T c712 = 240.0*c2;
  const T c713 = Bs*c707;
  const T c714 = 120.0*c28;
  const T c715 = 72.0*c602;
  const T c716 = 240.0*c696;
  const T c717 = 240.0*c532;
  const T c718 = 20.0*c365;
  const T c719 = PA1*c124;
  const T c720 = c48*c560;
  const T c721 = Bs*c573;
  const T c722 = 48.0*c28;
  const T c723 = 48.0*c372;
  const T c724 = 48.0*c435;
  const T c725 = 36.0*c395;
  const T c726 = 144.0*c2;
  const T c727 = 96.0*c28;
  const T c728 = 16.0*c365;
  const T c729 = 48.0*c438;
  const T c730 = 36.0*c602;
  const T c731 = c116*c577;
  const T c732 = c507 + c705;
  const T c733 = 12.0*c602;
  const T c734 = c580*d;
  const T c735 = c116*c150;
  const T c736 = c657 + c732;
  K[0] = to_dbl(c5);
  K[1] = to_dbl(c8);
  K[2] = to_dbl(-(T(1.0) / T(5040.0))*c4);
  K[3] = to_dbl(-(T(1.0) / T(720.0))*c4);
  K[4] = to_dbl(-(T(1.0) / T(120.0))*c4);
  K[5] = to_dbl(-(T(1.0) / T(24.0))*c4);
  K[6] = to_dbl(-d*(((T(1.0) / T(6.0)))*c1 + c2 + ((T(11.0) / T(6.0)))*p + 1.0));
  K[7] = to_dbl(-d*(((T(1.0) / T(2.0)))*c1 + c9 + ((T(11.0) / T(2.0)))*p + 3.0));
  K[8] = to_dbl(c5);
  K[9] = to_dbl(c5);
  K[10] = to_dbl(c8);
  K[11] = to_dbl(c13);
  K[12] = to_dbl(c19);
  K[13] = to_dbl(((T(1.0) / T(5040.0)))*c13);
  K[14] = to_dbl(((T(1.0) / T(720.0)))*c13);
  K[15] = to_dbl(((T(1.0) / T(120.0)))*c13);
  K[16] = to_dbl(((T(1.0) / T(24.0)))*c13);
  K[17] = to_dbl(((T(1.0) / T(6.0)))*c13);
  K[18] = to_dbl(c12*(((T(1.0) / T(2.0)))*c2 + ((T(3.0) / T(2.0)))*p + 1.0));
  K[19] = to_dbl(c13);
  K[20] = to_dbl(c13);
  K[21] = to_dbl(c19);
  K[22] = to_dbl(c22);
  K[23] = to_dbl(c22);
  K[24] = to_dbl(c26);
  K[25] = to_dbl(c26);
  K[26] = to_dbl(((T(1.0) / T(5040.0)))*c22);
  K[27] = to_dbl(-(T(1.0) / T(5040.0))*c21*(-c27 + c33 - 27.0*p - 16.0));
  K[28] = to_dbl(-(T(1.0) / T(720.0))*c21*(c33 - c34 - c35 - 14.0));
  K[29] = to_dbl(-(T(1.0) / T(120.0))*c21*(c33 - c36 - c37 - 12.0));
  K[30] = to_dbl(-(T(1.0) / T(24.0))*c21*(c33 - c38 - c39 - 10.0));
  K[31] = to_dbl(-(T(1.0) / T(6.0))*c21*(c33 - c40 + c41 - 8.0));
  K[32] = to_dbl(-(T(1.0) / T(2.0))*c21*(c33 - c42 - c43 - 6.0));
  K[33] = to_dbl(c22);
  K[34] = to_dbl(-c21*(c33 - c44 - c45 - 4.0));
  K[35] = to_dbl(c46);
  K[36] = to_dbl(c46);
  K[37] = to_dbl(c54);
  K[38] = to_dbl(c54);
  K[39] = to_dbl(-(T(1.0) / T(5040.0))*c53);
  K[40] = to_dbl(-(T(1.0) / T(720.0))*c53);
  K[41] = to_dbl(-(T(1.0) / T(120.0))*c53);
  K[42] = to_dbl(-(T(1.0) / T(24.0))*c53);
  K[43] = to_dbl(-(T(1.0) / T(6.0))*c53);
  K[44] = to_dbl(-(T(1.0) / T(2.0))*c53);
  K[45] = to_dbl(c54);
  K[46] = to_dbl(c54);
  K[47] = to_dbl(c54);
  K[48] = to_dbl(c57);
  K[49] = to_dbl(c57);
  K[50] = to_dbl(c72);
  K[51] = to_dbl(c72);
  K[52] = to_dbl(-(T(1.0) / T(5040.0))*c56);
  K[53] = to_dbl(c80*(-c73 + c75 + c79));
  K[54] = to_dbl(((T(1.0) / T(720.0)))*c71*(c79 + c81 + c83));
  K[55] = to_dbl(((T(1.0) / T(120.0)))*c71*(-7.0*c16 + c79 + c85));
  K[56] = to_dbl(((T(1.0) / T(24.0)))*c71*(-5.0*c58 + c79 + c87));
  K[57] = to_dbl(((T(1.0) / T(6.0)))*c71*(-5.0*c16 + c79 + c89));
  K[58] = to_dbl(((T(1.0) / T(2.0)))*c71*(-c23 + c79 + c91));
  K[59] = to_dbl(c57);
  K[60] = to_dbl(c71*(c79 + c93 - c94));
  K[61] = to_dbl(c95);
  K[62] = to_dbl(c95);
  K[63] = to_dbl(c97);
  K[64] = to_dbl(c97);
  K[65] = to_dbl(c104);
  K[66] = to_dbl(c104);
  K[67] = to_dbl(c122);
  K[68] = to_dbl(c122);
  K[69] = to_dbl(-(T(1.0) / T(5040.0))*c96);
  K[70] = to_dbl(c80*(-c0 + c106 + c123 + c15 + c31 + c92 - 9.0));
  K[71] = to_dbl(-c139*(c124 - c125 - c126 - 16.0*c127 - c128 - c129 - c130 + c138 + 56.0*p));
  K[72] = to_dbl(-c142*(c138 - c140 + c141 - 16.0*c17 + 56.0*c2 - 16.0*c28 + 42.0*p));
  K[73] = to_dbl(-c147*(-c12*c42 + c138 + c143 - c144 - c145 - c146 - 14.0*c17 + 42.0*c2 - 14.0*c28));
  K[74] = to_dbl(-c153*(-10.0*c101 - 10.0*c127 + c138 + c148 + c149 - c150 - c151 - c152 - 10.0*c58));
  K[75] = to_dbl(-c156*(-8.0*c127 + c138 + c154 - c155 - 10.0*c16 - 10.0*c17 + c42 - c74*d + c75));
  K[76] = to_dbl(c97);
  K[77] = to_dbl(-c160*(-c12*c6 + c138 + c157 - c158 - c159 - 8.0*c17 + c6 + c83 + c85));
  K[78] = to_dbl(c103*(c102 + c127 + c161));
  K[79] = to_dbl(-c121*(c100 + c138 - c162 - c163 - c164 - c165 + c3 + c87 + c89));
  K[80] = to_dbl(c168);
  K[81] = to_dbl(c168);
  K[82] = to_dbl(c173);
  K[83] = to_dbl(c173);
  K[84] = to_dbl(((T(1.0) / T(5040.0)))*c173);
  K[85] = to_dbl(((T(1.0) / T(720.0)))*c173);
  K[86] = to_dbl(((T(1.0) / T(120.0)))*c173);
  K[87] = to_dbl(((T(1.0) / T(24.0)))*c173);
  K[88] = to_dbl(((T(1.0) / T(6.0)))*c173);
  K[89] = to_dbl(((T(1.0) / T(2.0)))*c173);
  K[90] = to_dbl(c173);
  K[91] = to_dbl(c173);
  K[92] = to_dbl(c173);
  K[93] = to_dbl(c174);
  K[94] = to_dbl(c174);
  K[95] = to_dbl(c179);
  K[96] = to_dbl(c179);
  K[97] = to_dbl(((T(1.0) / T(5040.0)))*c174);
  K[98] = to_dbl(-c139*(c181 + c183 - 7.0*c47));
  K[99] = to_dbl(-c142*(c183 + c185 + c187));
  K[100] = to_dbl(-c147*(c183 + c189 + c191));
  K[101] = to_dbl(-c153*(c183 - c192 + c193));
  K[102] = to_dbl(-c156*(c183 + c195 + c197));
  K[103] = to_dbl(-c160*(c183 + c199 + c201));
  K[104] = to_dbl(c174);
  K[105] = to_dbl(-c121*(c183 + c202 - c47));
  K[106] = to_dbl(c203);
  K[107] = to_dbl(c203);
  K[108] = to_dbl(c204);
  K[109] = to_dbl(c204);
  K[110] = to_dbl(c209);
  K[111] = to_dbl(c209);
  K[112] = to_dbl(c222);
  K[113] = to_dbl(c222);
  K[114] = to_dbl(((T(1.0) / T(5040.0)))*c204);
  K[115] = to_dbl(-c139*(2.0*G1 + c117 + c224 + 2.0*c50));
  K[116] = to_dbl(((T(1.0) / T(5040.0)))*c221*(56.0*c16 + c227 - 14.0*c47 - 16.0*c49 - 14.0*c61 - 16.0*c62 - 14.0*c66 - 16.0*c77));
  K[117] = to_dbl(((T(1.0) / T(720.0)))*c221*(42.0*c16 + c227 + c228 + c231 - 14.0*c49 - 14.0*c62 - 14.0*c77));
  K[118] = to_dbl(((T(1.0) / T(120.0)))*c221*(c227 + c233 + c234 + c236 - 10.0*c47 - 10.0*c61 - 10.0*c66));
  K[119] = to_dbl(((T(1.0) / T(24.0)))*c221*(c227 - c237*d - c237 + c238 - 8.0*c47 - 8.0*c61 - 8.0*c66 - 10.0*c77));
  K[120] = to_dbl(((T(1.0) / T(6.0)))*c221*(c151 - c180*d + c181 + c227 + c240 + c242 - 8.0*c77));
  K[121] = to_dbl(c204);
  K[122] = to_dbl(((T(1.0) / T(2.0)))*c221*(c193 - c214 + c227 + c244 + c246 + c247 - 4.0*c61));
  K[123] = to_dbl(-c208*(c207 + c76));
  K[124] = to_dbl(c221*(-c110 - c114 - c194*d + c195 + c201 + c227 - c249 + c60));
  K[125] = to_dbl(c250);
  K[126] = to_dbl(c250);
  K[127] = to_dbl(c121);
  K[128] = to_dbl(c121);
  K[129] = to_dbl(c252);
  K[130] = to_dbl(c252);
  K[131] = to_dbl(c256);
  K[132] = to_dbl(c256);
  K[133] = to_dbl(c268);
  K[134] = to_dbl(c268);
  K[135] = to_dbl(c139);
  K[136] = to_dbl(-c139*(3.0*c12 + c269 + c270 - 10.0));
  K[137] = to_dbl(((T(1.0) / T(1680.0)))*c221*(-Bs*c44 - c12*c44 + c224 + c271 + c36));
  K[138] = to_dbl(-c282*(-G1*c36 - PA1*c36 - c116*c36 - 24.0*c134 + 168.0*c16 + 168.0*c17 - 336.0*c2 - c272 - 42.0*c273 - c274 - 48.0*c275 + 168.0*c28 + c281 - 42.0*c61 - 48.0*c62 - 48.0*c77));
  K[139] = to_dbl(-c286*(-G1*c38 + 126.0*c16 + 126.0*c17 - 210.0*c2 - 36.0*c273 + 126.0*c28 + c281 - c283 - c284 + c285 - 36.0*c61));
  K[140] = to_dbl(-c295*(-18.0*c113 - 18.0*c134 + 90.0*c17 + 90.0*c28 + c281 - c287 - c288 + c289 - c290 - c292 + c294 - 36.0*c62));
  K[141] = to_dbl(-c303*(-PA1*c300 - c116*c42 - 15.0*c134 + c231 - 30.0*c275 + c281 - c296 + c297 + c298 - c299 - c301 + c302 - 15.0*c49 - 30.0*c62 - 30.0*c77));
  K[142] = to_dbl(c121);
  K[143] = to_dbl(-c312*(-12.0*c113 - c116*c44 + c233 - c263 - 24.0*c275 + c281 - c283*d + c291 - c304 - c305 - c306 + c307 + c308 - c309 - c310 - c311));
  K[144] = to_dbl(-c161*c251);
  K[145] = to_dbl(-c316*(-Bs*c129 - c116*c6 + c125 + c126 + c129 - c131*c35 + c228 - c230*d + c242 + c281 - c313 - c314 - c315 + c41 - 18.0*c62));
  K[146] = to_dbl(c255*(c167 + c317));
  K[147] = to_dbl(-c267*(c162 + c164 + c197 + c236 + c240 + c247 + c281 - c318 - c319 - c320 - c321 - c322 - c323));
  K[148] = to_dbl(c324);
  K[149] = to_dbl(c324);
  K[150] = to_dbl(c332);
  K[151] = to_dbl(c339);
  K[152] = to_dbl(-(T(1.0) / T(5040.0))*c331);
  K[153] = to_dbl(-(T(1.0) / T(720.0))*c331);
  K[154] = to_dbl(-(T(1.0) / T(120.0))*c331);
  K[155] = to_dbl(-(T(1.0) / T(24.0))*c331);
  K[156] = to_dbl(-(T(1.0) / T(6.0))*c331);
  K[157] = to_dbl(-(T(1.0) / T(2.0))*c331);
  K[158] = to_dbl(c332);
  K[159] = to_dbl(c332);
  K[160] = to_dbl(c339);
  K[161] = to_dbl(c341);
  K[162] = to_dbl(c341);
  K[163] = to_dbl(c267*c350);
  K[164] = to_dbl(c266*c364);
  K[165] = to_dbl(-(T(1.0) / T(5040.0))*c340);
  K[166] = to_dbl(c282*(-5.0*c333 - 7.0*c365 + c366 + c369));
  K[167] = to_dbl(c286*(c330 + c338 + c348 + c368 + c371 - 15.0*c372));
  K[168] = to_dbl(c295*(c353 - 5.0*c365 + c369 + c374));
  K[169] = to_dbl(c303*(c369 + c376 + c377 + c379));
  K[170] = to_dbl(c312*(c369 + c380 + c382 + c384));
  K[171] = to_dbl(c316*(-c342 + c369 + c386));
  K[172] = to_dbl(c341);
  K[173] = to_dbl(c267*(c333 - c365 + c369));
  K[174] = to_dbl(c267*(c350 + c368));
  K[175] = to_dbl(c266*(Bs*c354 + c127*c175 + c364 + c387));
  K[176] = to_dbl(c389);
  K[177] = to_dbl(c389);
  K[178] = to_dbl(c392);
  K[179] = to_dbl(c392);
  K[180] = to_dbl(-c410*c412);
  K[181] = to_dbl(-c411*c426);
  K[182] = to_dbl(-(T(1.0) / T(5040.0))*c388);
  K[183] = to_dbl(c282*(c414*c66 + c428 + c429));
  K[184] = to_dbl(-(T(1.0) / T(5040.0))*c412*(G1*c149 - 10.0*c333 - 10.0*c395 - c430*d - c430 + c432 + c434 + c437 - 14.0*c438 - 10.0*c439 + c446 + c449));
  K[185] = to_dbl(-(T(1.0) / T(720.0))*c412*(G1*c154 - G1*c155 - G1*c234 - PA1*c238 - c234*c63 - 8.0*c439 + c446 + c451 + c452 + c455 + c458));
  K[186] = to_dbl(-(T(1.0) / T(120.0))*c412*(-G1*c158 - PA1*c140 + c421 + 20.0*c431 - 10.0*c438 + c446 - c459*d - c459 - 24.0*c460 + c462 + c465));
  K[187] = to_dbl(-(T(1.0) / T(24.0))*c412*(c363 + c407 - 8.0*c438 + c444 - c466*d + c469 + c470 + c472 + c473));
  K[188] = to_dbl(-(T(1.0) / T(6.0))*c412*(-G1*c166 - PA1*c82 + c409 + c419 + c444 - c474 - c475 + c476 + c479 + c480));
  K[189] = to_dbl(c389);
  K[190] = to_dbl(-(T(1.0) / T(2.0))*c412*(-c248*c431 - c375*d + c376 + c384 - c393 - c394 + 2.0*c431 + c446 - c482));
  K[191] = to_dbl(c391*(c182 + c390));
  K[192] = to_dbl(-c412*(-c30*c431 + c345 - c385*d + c386 + c397 + c418 + c445 + c474));
  K[193] = to_dbl(-c412*(c410 + c444 + c482));
  K[194] = to_dbl(-c411*(c111*c442 + c116*c245 + c123*c212 + c134*c63 + c248*c484 + c426 + c441*d + c443 + c481*c63 + c483));
  K[195] = to_dbl(c486);
  K[196] = to_dbl(c486);
  K[197] = to_dbl(c488);
  K[198] = to_dbl(c488);
  K[199] = to_dbl(c491);
  K[200] = to_dbl(c491);
  K[201] = to_dbl(c512*c514);
  K[202] = to_dbl(c513*c536);
  K[203] = to_dbl(-(T(1.0) / T(5040.0))*c485);
  K[204] = to_dbl(c282*(Bs*c23 + c537*d + c538 + 3.0*c63 + c81));
  K[205] = to_dbl(-(T(1.0) / T(1680.0))*c412*(c144 - c184*d + c226 - c416 + c429 - 5.0*c66 - 7.0*c77));
  K[206] = to_dbl(((T(1.0) / T(5040.0)))*c514*(G1*c540 - G2*c40 + PA1*c540 - PA2*c40 - c111*c539 - c143*c175 - c143*c544 - c212*c40 - c307*c63 - 54.0*c372 - 54.0*c435 - 210.0*c532 - c539 + c540*c63 + c541*d + c541 - c542*d - c542 - 42.0*c543 - 21.0*c545 + 126.0*c546 + c547 - 54.0*c548 + c556 + c557));
  K[207] = to_dbl(((T(1.0) / T(720.0)))*c514*(-G0*c560 + G1*c296 - G1*c563 + PA1*c296 - PA1*c302 - PA1*c563 - c111*c561 - c289*c63 + c296*c63 - 45.0*c372 - 45.0*c435 - 18.0*c545 + 90.0*c546 - 45.0*c548 + c556 + c559 - c561 + c562*d + c562 - c563*c63 - c565*d + c566 + c567 + c570));
  K[208] = to_dbl(((T(1.0) / T(120.0)))*c514*(Bs*c571 - G2*c44 - PA2*c44 - c111*c572 - c175*c38 - c38*c544 - 48.0*c484 + c509 + c526 - 30.0*c543 - 15.0*c545 + c555 - c571 - c572 + c573*d + c573 - c574*d - c574 + c575*c63 + c579 + c581 + c583));
  K[209] = to_dbl(((T(1.0) / T(24.0)))*c514*(-G1*c125 + G1*c516 - PA1*c125 + PA1*c516 - c125*c63 - 27.0*c372 - 27.0*c435 - 54.0*c460 + c462 - c492 + c516*c63 + 36.0*c546 - 27.0*c548 + c556 + c564*d + c587 + c590 + c592));
  K[210] = to_dbl(c486);
  K[211] = to_dbl(((T(1.0) / T(6.0)))*c514*(-PA1*c126 - c111*c593 + c353 - c355 + c421 + c470 - 24.0*c484 + c511 - c522 + c535 - 9.0*c545 + 18.0*c546 + c555 - c561*d - c593 + c594*d + c594 - c595*d - c595 - c596 - c597 + c599));
  K[212] = to_dbl(c487*(c206 + c76));
  K[213] = to_dbl(((T(1.0) / T(2.0)))*c514*(Bs*c533 - c111*c370 - c131*c533 + c379 - c454*d + c455 + c469 + c477*d + c480 - c493 - c494 - c495 + c556 - c600 - c601));
  K[214] = to_dbl(-c490*(c227 + c489));
  K[215] = to_dbl(c514*(-c111*c381 - c270*c602 + c352 + c355 + c382 + c420 - c478*d + c479 + c534 + c556 + c596));
  K[216] = to_dbl(c514*(c467 + c512 + c555 + c600 + c601));
  K[217] = to_dbl(c513*(Bs*c518 + Bs*c519 + 9.0*PA1*c275 + c111*c551 + c116*c162 + c116*c277 + c116*c314 + c12*c603 + c159*c212 + c159*c417 + c175*c29 + c235*c605 + c258*c549 + c259*c29 + c269*c543 + c29*c544 + c536 + c550*d + c554 + c600*d + c601 + 6.0*c604));
  K[218] = to_dbl(c606);
  K[219] = to_dbl(c606);
  K[220] = to_dbl(c608);
  K[221] = to_dbl(c608);
  K[222] = to_dbl(c610);
  K[223] = to_dbl(c610);
  K[224] = to_dbl(c613);
  K[225] = to_dbl(c613);
  K[226] = to_dbl(-c662*c664);
  K[227] = to_dbl(-c663*c676);
  K[228] = to_dbl(-c282);
  K[229] = to_dbl(c282*(c248 + c414 + c607 - 7.0));
  K[230] = to_dbl(-(T(1.0) / T(2520.0))*c412*(Bs*c84 + PA1*c552 + c10*c131 + c111*c537 + 3.0*c116 + c12*c159 + c141 + c279*c6 + c36 + c538 + 6.0*c63));
  K[231] = to_dbl(((T(1.0) / T(2520.0)))*c514*(G0*c284 + G0*c311 + G1*c159 + PA1*c159 + PA2*c552 + c108*c276 + c111*c677 + c111*c678 + c116*c159 + c12*c283 + c12*c318 + c131*c86 + c145*c278 + c145*c63 + c146*c605 + c150*c279 + 84.0*c16 + 84.0*c17 + 6.0*c175 - 105.0*c2 + c210*c258 + 6.0*c212 + c258*c427 + 2.0*c259 + c262 + 84.0*c28 + c285 + c294 + 6.0*c417 + c428 + c677 + c678*d));
  K[232] = to_dbl(-(T(1.0) / T(5040.0))*c664*(-Bs*c679 - G0*c679 + G1*c681 - G1*c685 - G2*c148 + PA1*c681 - PA1*c685 - PA2*c148 + c111*c682 - c111*c683 - c111*c689 + c116*c681 - c116*c685 - c12*c679 - c148*c259 - c258*c680 - c258*c686 + c278*c690 - c278*c693 - 60.0*c439 - 216.0*c460 - 216.0*c484 - 84.0*c543 - 84.0*c545 + 504.0*c546 - 216.0*c548 + 252.0*c602 - 28.0*c603 + c63*c690 - c63*c693 + c679 - c680 + c682 - c683 - 84.0*c684 - c686 - c687 - 60.0*c688 - c689*d + 504.0*c691 - 168.0*c692 - 144.0*c694 - 84.0*c695 + 504.0*c697 + c706));
  K[233] = to_dbl(-(T(1.0) / T(720.0))*c664*(G0*c713 + G1*c287 - G1*c297 - G1*c298 + PA1*c287 - PA1*c297 - PA1*c298 + c111*c708 - c111*c710 - c111*c711 + c116*c287 - c116*c297 + c12*c713 - c12*c715 - c258*c588 - c276*c584 + c278*c712 - c278*c714 + c279*c707 - c298*c523 - 16.0*c333 - 16.0*c395 - 48.0*c439 - 180.0*c460 - 180.0*c484 - 480.0*c517 - 480.0*c532 - 72.0*c543 - 72.0*c545 - 180.0*c548 - c588 + 180.0*c602 - 16.0*c614 + c63*c712 - c63*c714 - 72.0*c684 - 48.0*c688 - 144.0*c692 - 120.0*c694 - 480.0*c696 + c706 + c707 + c708 - c710 - c711*d));
  K[234] = to_dbl(-(T(1.0) / T(120.0))*c664*(Bs*c717 - Bs*c720*d + G1*c124 - G1*c722 - PA1*c722 + c111*c720 - c111*c721 - c111*c725 + c116*c124 + c12*c716 + c131*c560 - c131*c571 - c258*c718 - c258*c723 + c278*c726 - c278*c727 - c28*c666 - 36.0*c439 - 144.0*c460 - 144.0*c484 - 240.0*c517 - 60.0*c543 - 144.0*c548 + c560 + c567 - 20.0*c603 + c63*c726 - c63*c727 + c667 - 60.0*c684 - 36.0*c688 + 240.0*c691 - 96.0*c694 - 60.0*c695 + c706 - c716 - c717 - c718 + c719 + c720 - c721 - c723 - c724 - c725*d));
  K[235] = to_dbl(c606);
  K[236] = to_dbl(-(T(1.0) / T(24.0))*c664*(-c111*c687 + c111*c709 - c111*c729 - c116*c307 + c116*c575 + c124*c63 - c258*c728 - c278*c685 + c458 - 108.0*c484 - 96.0*c517 - 96.0*c532 - 48.0*c543 - 48.0*c545 + 144.0*c546 - 108.0*c548 + c557 + c570 + c583 + c584 - 16.0*c603 - c615 - c616 - c617 - c63*c685 - c672 - 48.0*c684 + 144.0*c691 - 96.0*c692 - 48.0*c695 - 96.0*c696 + 144.0*c697 + c706 + c709 + c715 + c719*d - c728 - c729));
  K[237] = to_dbl(4.0*c267*(Bs + c98));
  K[238] = to_dbl(-(T(1.0) / T(6.0))*c664*(c111*c564 - c111*c565 - c111*c724 - c12*c730 - c258*c450 - c258*c464 - c278*c722 + c437 + c465 + c473 + c531 - 36.0*c545 + 72.0*c546 + c547 - 72.0*c548 + c566 + c581 - c585 + c592 - 12.0*c603 + c623 + c624 + c625 + c626 + c627 + c628 + c629 - c63*c722 + c630 + c631 + c632 + c634 + c635 + c637 + c639 + c640 + c641 + c642 + c643 + c644 + c645 + c646 + c647 + c648 + c649 + c650 + c651 + c653 + c654 + c656 + c661 + c675 - 36.0*c684 + 72.0*c697 + c709*d - c710*d + c730 - c731 + c732));
  K[239] = to_dbl(-c317*c609);
  K[240] = to_dbl(-(T(1.0) / T(2.0))*c664*(Bs*c591 + c111*c453 - c111*c586 - c111*c588 + c12*c585 - c131*c591 - c131*c673 - c258*c466 + c279*c584 + c472 + c476 + c579 + c590 + c599 - 8.0*c603 - c618 - c619 - c729*d + c733 - c734 - c735 + c736));
  K[241] = to_dbl(c612*(c280 + c611));
  K[242] = to_dbl(-c664*(-c111*c454 - c12*c733 - c258*c375 + c376 + c399 + c455 + c499 + c527 + c558 + c587 - 4.0*c603 + c674 + c736));
  K[243] = to_dbl(-c664*(c447 + c578 + c598 + c662 + c704 + c734 + c735));
  K[244] = to_dbl(-c663*(12.0*c111*c440 + c111*c580 + c111*c701 + c116*c305 + c145*c259 + c145*c520 + c163*c175 + c163*c544 + c235*c63 + c235*c670 + c248*c543 + 12.0*c258*c367 + c258*c447 + c258*c702 + c278*c318 + c290*c523 + c290*c605 + c305*c605 + c309*c521 + c318*c523 + c318*c63 + c440*c552 + c528*c699 + c528*c700 + c553 + c597 + c598*d + c676 + c698*d + c703 + c731));
}

// order 3, the unique components in the order of the kernel (111, 112, 113,
// 122, 123, 133, 222, 223, 233, 333) at Z, pC, A0 and im = 1/mean
inline void gengamma2_eval3(const double* K, double Z, double pC, double A0,
                            double im, double* u) {
    if (std::fabs(Z) < 1.0) {
      const double R = gg2_R(Z);
      const double u0 = std::pow(Z, 7);
      const double u1 = std::pow(Z, 6);
      const double u2 = std::pow(Z, 5);
      const double u3 = std::pow(Z, 4);
      const double u4 = std::pow(Z, 3);
      const double u5 = std::pow(Z, 2);
      const double u6 = R*Z;
      const double u7 = std::pow(Z, 8);
      const double u8 = std::pow(Z, 9);
      const double u9 = R*u5;
      u[0] = (K[10] + K[2]*u0 + K[3]*u1 + K[4]*u2 + K[5]*u3 + K[6]*u4 + K[7]*u5 + K[8]*Z + K[9]*R) * std::pow(im, 3);
      u[1] = (K[13]*u0 + K[14]*u1 + K[15]*u2 + K[16]*u3 + K[17]*u4 + K[18]*u5 + K[19]*Z + K[20]*R + K[21]) * std::pow(im, 2);
      u[2] = (K[26]*u7 + K[27]*u0 + K[28]*u1 + K[29]*u2 + K[30]*u3 + K[31]*u4 + K[32]*u5 + K[33]*u6 + K[34]*Z + K[35]*R + K[36]) * std::pow(im, 2);
      u[3] = (K[39]*u0 + K[40]*u1 + K[41]*u2 + K[42]*u3 + K[43]*u4 + K[44]*u5 + K[45]*Z + K[46]*R + K[47]) * im;
      u[4] = (K[52]*u7 + K[53]*u0 + K[54]*u1 + K[55]*u2 + K[56]*u3 + K[57]*u4 + K[58]*u5 + K[59]*u6 + K[60]*Z + K[61]*R + K[62]) * im;
      u[5] = (K[69]*u8 + K[70]*u7 + K[71]*u0 + K[72]*u1 + K[73]*u2 + K[74]*u3 + K[75]*u4 + K[76]*u9 + K[77]*u5 + K[78]*u6 + K[79]*Z + K[80]*R + K[81]) * im;
      u[6] = (K[84]*u0 + K[85]*u1 + K[86]*u2 + K[87]*u3 + K[88]*u4 + K[89]*u5 + K[90]*Z + K[91]*R + K[92]);
      u[7] = (K[100]*u2 + K[101]*u3 + K[102]*u4 + K[103]*u5 + K[104]*u6 + K[105]*Z + K[106]*R + K[107] + K[97]*u7 + K[98]*u0 + K[99]*u1);
      u[8] = (K[114]*u8 + K[115]*u7 + K[116]*u0 + K[117]*u1 + K[118]*u2 + K[119]*u3 + K[120]*u4 + K[121]*u9 + K[122]*u5 + K[123]*u6 + K[124]*Z + K[125]*R + K[126]);
      u[9] = (K[135]*std::pow(Z, 10) + K[136]*u8 + K[137]*u7 + K[138]*u0 + K[139]*u1 + K[140]*u2 + K[141]*u3 + K[142]*R*u4 + K[143]*u4 + K[144]*u9 + K[145]*u5 + K[146]*u6 + K[147]*Z + K[148]*R + K[149]);
    } else {
      const double X = pC + A0;
      const double Q = std::expm1(Z);
      const double u0 = Q*X;
      const double u1 = std::pow(X, 2);
      const double u2 = Q*u1;
      const double u3 = std::pow(X, 3);
      u[0] = (K[0]*Q + K[1]) * std::pow(im, 3);
      u[1] = (K[11]*Q + K[12]) * std::pow(im, 2);
      u[2] = (K[22]*u0 + K[23]*X + K[24]*Q + K[25]) * std::pow(im, 2);
      u[3] = (K[37]*Q + K[38]) * im;
      u[4] = (K[48]*u0 + K[49]*X + K[50]*Q + K[51]) * im;
      u[5] = (K[63]*u2 + K[64]*u1 + K[65]*u0 + K[66]*X + K[67]*Q + K[68]) * im;
      u[6] = (K[82]*Q + K[83]);
      u[7] = (K[93]*u0 + K[94]*X + K[95]*Q + K[96]);
      u[8] = (K[108]*u2 + K[109]*u1 + K[110]*u0 + K[111]*X + K[112]*Q + K[113]);
      u[9] = (K[127]*Q*u3 + K[128]*u3 + K[129]*u2 + K[130]*u1 + K[131]*u0 + K[132]*X + K[133]*Q + K[134]);
    }
}

// order 4, the unique components in the order of the kernel (1111, 1112,
// 1113, 1122, 1123, 1133, 1222, 1223, 1233, 1333, 2222, 2223, 2233, 2333,
// 3333)
inline void gengamma2_eval4(const double* K, double Z, double pC, double A0,
                            double im, double* u) {
    if (std::fabs(Z) < 1.0) {
      const double R = gg2_R(Z);
      const double u0 = std::pow(Z, 7);
      const double u1 = std::pow(Z, 6);
      const double u2 = std::pow(Z, 5);
      const double u3 = std::pow(Z, 4);
      const double u4 = std::pow(Z, 3);
      const double u5 = std::pow(Z, 2);
      const double u6 = R*Z;
      const double u7 = std::pow(Z, 8);
      const double u8 = std::pow(Z, 9);
      const double u9 = R*u5;
      const double u10 = std::pow(Z, 10);
      const double u11 = R*u4;
      u[0] = (K[10] + K[2]*u0 + K[3]*u1 + K[4]*u2 + K[5]*u3 + K[6]*u4 + K[7]*u5 + K[8]*Z + K[9]*R) * std::pow(im, 4);
      u[1] = (K[13]*u0 + K[14]*u1 + K[15]*u2 + K[16]*u3 + K[17]*u4 + K[18]*u5 + K[19]*Z + K[20]*R + K[21]) * std::pow(im, 3);
      u[2] = (K[26]*u7 + K[27]*u0 + K[28]*u1 + K[29]*u2 + K[30]*u3 + K[31]*u4 + K[32]*u5 + K[33]*u6 + K[34]*Z + K[35]*R + K[36]) * std::pow(im, 3);
      u[3] = (K[39]*u0 + K[40]*u1 + K[41]*u2 + K[42]*u3 + K[43]*u4 + K[44]*u5 + K[45]*Z + K[46]*R + K[47]) * std::pow(im, 2);
      u[4] = (K[52]*u7 + K[53]*u0 + K[54]*u1 + K[55]*u2 + K[56]*u3 + K[57]*u4 + K[58]*u5 + K[59]*u6 + K[60]*Z + K[61]*R + K[62]) * std::pow(im, 2);
      u[5] = (K[69]*u8 + K[70]*u7 + K[71]*u0 + K[72]*u1 + K[73]*u2 + K[74]*u3 + K[75]*u4 + K[76]*u9 + K[77]*u5 + K[78]*u6 + K[79]*Z + K[80]*R + K[81]) * std::pow(im, 2);
      u[6] = (K[84]*u0 + K[85]*u1 + K[86]*u2 + K[87]*u3 + K[88]*u4 + K[89]*u5 + K[90]*Z + K[91]*R + K[92]) * im;
      u[7] = (K[100]*u2 + K[101]*u3 + K[102]*u4 + K[103]*u5 + K[104]*u6 + K[105]*Z + K[106]*R + K[107] + K[97]*u7 + K[98]*u0 + K[99]*u1) * im;
      u[8] = (K[114]*u8 + K[115]*u7 + K[116]*u0 + K[117]*u1 + K[118]*u2 + K[119]*u3 + K[120]*u4 + K[121]*u9 + K[122]*u5 + K[123]*u6 + K[124]*Z + K[125]*R + K[126]) * im;
      u[9] = (K[135]*u10 + K[136]*u8 + K[137]*u7 + K[138]*u0 + K[139]*u1 + K[140]*u2 + K[141]*u3 + K[142]*u11 + K[143]*u4 + K[144]*u9 + K[145]*u5 + K[146]*u6 + K[147]*Z + K[148]*R + K[149]) * im;
      u[10] = (K[152]*u0 + K[153]*u1 + K[154]*u2 + K[155]*u3 + K[156]*u4 + K[157]*u5 + K[158]*Z + K[159]*R + K[160]);
      u[11] = (K[165]*u7 + K[166]*u0 + K[167]*u1 + K[168]*u2 + K[169]*u3 + K[170]*u4 + K[171]*u5 + K[172]*u6 + K[173]*Z + K[174]*R + K[175]);
      u[12] = (K[182]*u8 + K[183]*u7 + K[184]*u0 + K[185]*u1 + K[186]*u2 + K[187]*u3 + K[188]*u4 + K[189]*u9 + K[190]*u5 + K[191]*u6 + K[192]*Z + K[193]*R + K[194]);
      u[13] = (K[203]*u10 + K[204]*u8 + K[205]*u7 + K[206]*u0 + K[207]*u1 + K[208]*u2 + K[209]*u3 + K[210]*u11 + K[211]*u4 + K[212]*u9 + K[213]*u5 + K[214]*u6 + K[215]*Z + K[216]*R + K[217]);
      u[14] = (K[228]*std::pow(Z, 11) + K[229]*u10 + K[230]*u8 + K[231]*u7 + K[232]*u0 + K[233]*u1 + K[234]*u2 + K[235]*R*u3 + K[236]*u3 + K[237]*u11 + K[238]*u4 + K[239]*u9 + K[240]*u5 + K[241]*u6 + K[242]*Z + K[243]*R + K[244]);
    } else {
      const double X = pC + A0;
      const double Q = std::expm1(Z);
      const double u0 = Q*X;
      const double u1 = std::pow(X, 2);
      const double u2 = Q*u1;
      const double u3 = std::pow(X, 3);
      const double u4 = Q*u3;
      const double u5 = std::pow(X, 4);
      u[0] = (K[0]*Q + K[1]) * std::pow(im, 4);
      u[1] = (K[11]*Q + K[12]) * std::pow(im, 3);
      u[2] = (K[22]*u0 + K[23]*X + K[24]*Q + K[25]) * std::pow(im, 3);
      u[3] = (K[37]*Q + K[38]) * std::pow(im, 2);
      u[4] = (K[48]*u0 + K[49]*X + K[50]*Q + K[51]) * std::pow(im, 2);
      u[5] = (K[63]*u2 + K[64]*u1 + K[65]*u0 + K[66]*X + K[67]*Q + K[68]) * std::pow(im, 2);
      u[6] = (K[82]*Q + K[83]) * im;
      u[7] = (K[93]*u0 + K[94]*X + K[95]*Q + K[96]) * im;
      u[8] = (K[108]*u2 + K[109]*u1 + K[110]*u0 + K[111]*X + K[112]*Q + K[113]) * im;
      u[9] = (K[127]*u4 + K[128]*u3 + K[129]*u2 + K[130]*u1 + K[131]*u0 + K[132]*X + K[133]*Q + K[134]) * im;
      u[10] = (K[150]*Q + K[151]);
      u[11] = (K[161]*u0 + K[162]*X + K[163]*Q + K[164]);
      u[12] = (K[176]*u2 + K[177]*u1 + K[178]*u0 + K[179]*X + K[180]*Q + K[181]);
      u[13] = (K[195]*u4 + K[196]*u3 + K[197]*u2 + K[198]*u1 + K[199]*u0 + K[200]*X + K[201]*Q + K[202]);
      u[14] = (K[218]*Q*u5 + K[219]*u5 + K[220]*u4 + K[221]*u3 + K[222]*u2 + K[223]*u1 + K[224]*u0 + K[225]*X + K[226]*Q + K[227]);
    }
}

// The buffer of the parts route, P = 3 in (mean, d, p). Each double-double
// quantity takes two slots (hi, lo; lo = 0 in the double branch). Index,
// quantity, lowest level:
//   [0] level   [1] mean   [2] d   [3] p
//   [4] 1 if the coefficients are formed in double (k0 < GG2_KC and
//       d < GG2_DC), 0 if in double-double
//   [5] pC = p log(y/mean) (1)   [6] A0 (1)   [7] Z0 = A0 + B0 (1)
//   [8,9]   Bs = psi(k0) - log(k0) in the branch's arithmetic (1)
//   [10,11] Gv0 = psi(k1) - psi(k0) (1)
//   [12,13] Gv1 = psi'(k1) - psi'(k0) (2)
//   [14,15] Gv2 (3)    [16,17] Gv3 (4)
//   [18,19] PAv1 = psi'(k0) (2)   [20,21] PAv2 (3)   [22,23] PAv3 (4)
// with k0 = d/p and k1 = (d + 1)/p. The expected parts leave pC, A0 and Z0
// as NaN and store Bs = 0, the Bs the expected coefficients are formed with.
inline void gg2p_put(double* b, int i, double v) { b[i] = v; b[i + 1] = 0.0; }
inline void gg2p_put(double* b, int i, const dd& v) { b[i] = v.hi; b[i + 1] = v.lo; }
inline void gg2p_get(const double* b, int i, double& v) { v = b[i]; }
inline void gg2p_get(const double* b, int i, dd& v) { v = dd(b[i], b[i + 1]); }

inline void gg2p_nan(double* b, int i) { b[i] = R_NaN; b[i + 1] = R_NaN; }

// the polygammas of the buffer, formed exactly as the kernels form them
inline void gengamma2_parts_poly(double dv, double pv, int level, double Bs,
                                 bool expected, double* b) {
  const double k0d = dv / pv;
  const bool dbl = k0d < GG2_KC && dv < GG2_DC;
  b[4] = dbl ? 1.0 : 0.0;
  for (int i = 8; i <= 22; i += 2) gg2p_nan(b, i);
  if (dbl) {
    const double k1d = (dv + 1.0) / pv;
    if (level >= 1) {
      gg2p_put(b, 8, expected ? 0.0 : Bs);
      gg2p_put(b, 10, R::psigamma(k1d, 0.0) - R::psigamma(k0d, 0.0));
    }
    if (level >= 2) {
      const double pa = R::psigamma(k0d, 1.0);
      gg2p_put(b, 18, pa);
      gg2p_put(b, 12, R::psigamma(k1d, 1.0) - pa);
    }
    if (level >= 3) {
      const double pa = R::psigamma(k0d, 2.0);
      gg2p_put(b, 20, pa);
      gg2p_put(b, 14, R::psigamma(k1d, 2.0) - pa);
    }
    if (level >= 4) {
      const double pa = R::psigamma(k0d, 3.0);
      gg2p_put(b, 22, pa);
      gg2p_put(b, 16, R::psigamma(k1d, 3.0) - pa);
    }
  } else {
    const dd d(dv), p(pv);
    const dd k0 = d / p, k1 = (d + 1.0) / p;
    if (level >= 1) {
      gg2p_put(b, 8, expected ? dd(0.0) : dd_B0(k0));
      gg2p_put(b, 10, dd_psi_diff(k0, k1, 0));
    }
    if (level >= 2) {
      gg2p_put(b, 18, dd_psi_n(k0, 1));
      gg2p_put(b, 12, dd_psi_diff(k0, k1, 1));
    }
    if (level >= 3) {
      gg2p_put(b, 20, dd_psi_n(k0, 2));
      gg2p_put(b, 14, dd_psi_diff(k0, k1, 2));
    }
    if (level >= 4) {
      gg2p_put(b, 22, dd_psi_n(k0, 3));
      gg2p_put(b, 16, dd_psi_diff(k0, k1, 3));
    }
  }
}

inline void gengamma2_parts(double y, const double* th, int level, double* b) {
  const double mv = th[0], dv = th[1], pv = th[2];
  b[0] = level; b[1] = mv; b[2] = dv; b[3] = pv;
  b[5] = b[6] = b[7] = R_NaN;
  double B0d = R_NaN;
  if (level >= 1) {
    const double A0 = gg2_A0(dv, pv);
    B0d = gg2_B0(dv / pv);
    const double dy = y - mv;
    b[5] = pv * (std::fabs(dy) < 0.5 * mv ? std::log1p(dy / mv) : std::log(y / mv));
    b[6] = A0;
    b[7] = A0 + B0d;
  }
  gengamma2_parts_poly(dv, pv, level, B0d, false, b);
}

inline void gengamma2_parts_expected(const double* th, int level, double* b) {
  const double mv = th[0], dv = th[1], pv = th[2];
  b[0] = level; b[1] = mv; b[2] = dv; b[3] = pv;
  b[5] = b[6] = b[7] = R_NaN;
  gengamma2_parts_poly(dv, pv, level, 0.0, true, b);
}

// the arrays the coefficient generators read, from the buffer: only the
// slots the kernel of that order sets (nG polygamma differences from Gv0, and
// the derivatives of psi from PAv1 up to nPA)
template <class T>
inline void gengamma2_parts_load(const double* b, int nG, int nPA, T& d, T& p,
                                 T* PAv, T* Gv, T& Bs) {
  d = T(b[2]); p = T(b[3]);
  for (int j = 0; j < 12; ++j) { PAv[j] = T(0.0); Gv[j] = T(0.0); }
  for (int j = 0; j < nG; ++j) gg2p_get(b, 10 + 2 * j, Gv[j]);
  for (int j = 1; j <= nPA; ++j) gg2p_get(b, 16 + 2 * j, PAv[j]);
  gg2p_get(b, 8, Bs);
}

inline GG2Data gengamma2_parts_data(const double* b, double& im) {
  im = 1.0 / b[1];
  const double pC = b[5];
  const double Z = pC + b[7];
  return gg2_data(Z, pC, b[6]);
}

inline void gengamma2_parts_grad(const double* b, double* out) {
  double im = 0.0;
  const GG2Data D = gengamma2_parts_data(b, im);
  double K[35] = {0};
  if (b[4] != 0.0) {
    double d, p, Bs, PAv[12], Gv[12];
    gengamma2_parts_load<double>(b, 1, 0, d, p, PAv, Gv, Bs);
    gengamma2_coef_score_mean<double>(d, p, PAv, Gv, Bs, K);
    gengamma2_coef_score_d<double>(d, p, PAv, Gv, Bs, K);
    gengamma2_coef_score_p<double>(d, p, PAv, Gv, Bs, K);
  } else {
    dd d, p, Bs, PAv[12], Gv[12];
    gengamma2_parts_load<dd>(b, 1, 0, d, p, PAv, Gv, Bs);
    gengamma2_coef_score_mean<dd>(d, p, PAv, Gv, Bs, K);
    gengamma2_coef_score_d<dd>(d, p, PAv, Gv, Bs, K);
    gengamma2_coef_score_p<dd>(d, p, PAv, Gv, Bs, K);
  }
  out[0] = gengamma2_score_mean(K, D, im);
  out[1] = gengamma2_score_d(K, D, im);
  out[2] = gengamma2_score_p(K, D, im);
}

inline void gengamma2_parts_hess(const double* b, double* out) {
  double im = 0.0;
  const GG2Data D = gengamma2_parts_data(b, im);
  double K[82] = {0};
  if (b[4] != 0.0) {
    double d, p, Bs, PAv[12], Gv[12];
    gengamma2_parts_load<double>(b, 2, 1, d, p, PAv, Gv, Bs);
    gengamma2_coef_hess_mean_mean<double>(d, p, PAv, Gv, Bs, K);
    gengamma2_coef_hess_d_d<double>(d, p, PAv, Gv, Bs, K);
    gengamma2_coef_hess_p_p<double>(d, p, PAv, Gv, Bs, K);
    gengamma2_coef_hessian_rest<double>(d, p, PAv, Gv, Bs, K);
  } else {
    dd d, p, Bs, PAv[12], Gv[12];
    gengamma2_parts_load<dd>(b, 2, 1, d, p, PAv, Gv, Bs);
    gengamma2_coef_hess_mean_mean<dd>(d, p, PAv, Gv, Bs, K);
    gengamma2_coef_hess_d_d<dd>(d, p, PAv, Gv, Bs, K);
    gengamma2_coef_hess_p_p<dd>(d, p, PAv, Gv, Bs, K);
    gengamma2_coef_hessian_rest<dd>(d, p, PAv, Gv, Bs, K);
  }
  out[0] = gengamma2_hess_mean_mean(K, D, im);
  out[4] = gengamma2_hess_d_d(K, D, im);
  out[8] = gengamma2_hess_p_p(K, D, im);
  out[1] = out[3] = gengamma2_hess_mean_d(K, D, im);
  out[2] = out[6] = gengamma2_hess_mean_p(K, D, im);
  out[5] = out[7] = gengamma2_hess_d_p(K, D, im);
}

// from the expected parts
inline void gengamma2_parts_info(const double* b, double* out) {
  const double im = 1.0 / b[1];
  double K[6] = {0};
  if (b[4] != 0.0) {
    double d, p, Bs, PAv[12], Gv[12];
    gengamma2_parts_load<double>(b, 1, 1, d, p, PAv, Gv, Bs);
    gengamma2_coef_expected_mean_mean<double>(d, p, PAv, Gv, Bs, K);
    gengamma2_coef_expected_d_d<double>(d, p, PAv, Gv, Bs, K);
    gengamma2_coef_expected_p_p<double>(d, p, PAv, Gv, Bs, K);
    gengamma2_coef_expected_hessian_rest<double>(d, p, PAv, Gv, Bs, K);
  } else {
    dd d, p, Bs, PAv[12], Gv[12];
    gengamma2_parts_load<dd>(b, 1, 1, d, p, PAv, Gv, Bs);
    gengamma2_coef_expected_mean_mean<dd>(d, p, PAv, Gv, Bs, K);
    gengamma2_coef_expected_d_d<dd>(d, p, PAv, Gv, Bs, K);
    gengamma2_coef_expected_p_p<dd>(d, p, PAv, Gv, Bs, K);
    gengamma2_coef_expected_hessian_rest<dd>(d, p, PAv, Gv, Bs, K);
  }
  out[0] = gengamma2_expected_mean_mean(K, im);
  out[4] = gengamma2_expected_d_d(K, im);
  out[8] = gengamma2_expected_p_p(K, im);
  out[1] = out[3] = gengamma2_expected_mean_d(K, im);
  out[2] = out[6] = gengamma2_expected_mean_p(K, im);
  out[5] = out[7] = gengamma2_expected_d_p(K);
}

inline void gengamma2_parts_d3(const double* b, double* out) {
  double K[150] = {0};
  if (b[4] != 0.0) {
    double d, p, Bs, PAv[12], Gv[12];
    gengamma2_parts_load<double>(b, 3, 2, d, p, PAv, Gv, Bs);
    gengamma2_coef_deriv3<double>(d, p, PAv, Gv, Bs, K);
  } else {
    dd d, p, Bs, PAv[12], Gv[12];
    gengamma2_parts_load<dd>(b, 3, 2, d, p, PAv, Gv, Bs);
    gengamma2_coef_deriv3<dd>(d, p, PAv, Gv, Bs, K);
  }
  const double im = 1.0 / b[1];
  double u[10];
  gengamma2_eval3(K, b[5] + b[7], b[5], b[6], im, u);
  sym_fill(3, 3, u, out);
}

inline void gengamma2_parts_d4(const double* b, double* out) {
  double K[245] = {0};
  if (b[4] != 0.0) {
    double d, p, Bs, PAv[12], Gv[12];
    gengamma2_parts_load<double>(b, 4, 3, d, p, PAv, Gv, Bs);
    gengamma2_coef_deriv4<double>(d, p, PAv, Gv, Bs, K);
  } else {
    dd d, p, Bs, PAv[12], Gv[12];
    gengamma2_parts_load<dd>(b, 4, 3, d, p, PAv, Gv, Bs);
    gengamma2_coef_deriv4<dd>(d, p, PAv, Gv, Bs, K);
  }
  const double im = 1.0 / b[1];
  double u[15];
  gengamma2_eval4(K, b[5] + b[7], b[5], b[6], im, u);
  sym_fill(3, 4, u, out);
}

} // namespace d7

#endif
