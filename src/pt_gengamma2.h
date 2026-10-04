#ifndef D7_PT_GENGAMMA2_H
#define D7_PT_GENGAMMA2_H

#include <Rcpp.h>
#include <cmath>

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

} // namespace d7

#endif
