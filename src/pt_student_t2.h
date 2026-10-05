#ifndef D7_PT_STUDENT_T2_H
#define D7_PT_STUDENT_T2_H

#include <Rcpp.h>
#include <cmath>

// The Student t in its location mu, standard deviation s and degrees of freedom nu:
// one function per component and order for the quantities the scalar
// registry reads, lifted from the generated kernels in student_t2.cpp (which
// records the derivation) and called by them and by the registry
// (d7_ccallable.cpp). With ik = 1/(nu - 2) and z = (y - mu)/s, the caller
// passes z, ik, D(q) for q = z^2 ik, and the quantities of nu alone
// t2_V*(nu), each of which carries a polygamma.

// D(q) = q/(1+q) - log1p(q). The direct form cancels as q -> 0, where D is
// -q^2/2. With w = q/(2+q), q/(1+q) = 2w/(1+w) and log1p(q) = 2 atanh(w), so
//   D = -2w^2/(1+w) - 2 sum_{j>=1} w^(2j+1)/(2j+1),
// a sum of terms of one sign; below q = 0.5 (w <= 0.2) twelve terms reach
// the last bit, and above it the direct form loses at most a few.
static inline double t2_D(double q) {
  if (q < 0.5) {
    const double w = q / (2.0 + q), w2 = w * w;
    double acc = 1.0 / 25.0;
    for (int j = 11; j >= 1; j--) acc = 1.0 / (2.0 * j + 1.0) + w2 * acc;
    return -2.0 * w2 / (1.0 + w) - 2.0 * w * w2 * acc;
  }
  return q / (1.0 + q) - std::log1p(q);
}

static inline double t2_V0(double nu) {
  if (nu >= 20.0) {
    const double h = 1.0 / nu;
    return std::pow(h, 2) * (-0.75 + h * (-2.0 + h * (-4.125 + h * (-8.0 + h * (-15.75 + h * (-32.0 + h * (-65.0625 + h * (-128.0 + h * (-248.25 + h * (-512.0 + h * (-1110.375 + h * (-2048.0 + h * (-2730.75 + h * (-8192.0 + h * (-45433.03125 + h * (-32768.0 + h * (735036.75 + h * (-131072.0 + h * (-28003466.625 + h * (-524288.0 + h * (1179480554.25 + h * (-2097152.0 + h * (-60528174355.6875 + h * (-8388608.0 + h * (3679400001321.75 + h * (-33554432.0 + h * (-261707677015447.88 + h * (-134217728.0 + h * (2.153141787236484e+16 + h * (-536870912.0 + h * (-2.0288775585910433e+18 + h * (-2147483648.0 + h * (2.1708009902194275e+20 + h * (-8589934592.0 + h * (-2.6173826968472997e+22 + h * (-34359738368.0 + h * (3.532414887686319e+24 + h * (-137438953472.0 + h * (-5.304203340686493e+26 + h * (-549755813888.0))))))))))))))))))))))))))))))))))))))));
  }
  const double PA0 = R::psigamma(0.5 * (nu + 1.0), 0.0);
  const double PB0 = R::psigamma(0.5 * nu, 0.0);
  return (1.0/2.0)*PA0 - 1.0/2.0*PB0 - (1.0/2.0)/(nu - 2);
}

static inline double t2_V1(double nu) {
  if (nu >= 20.0) {
    const double h = 1.0 / nu;
    return std::pow(h, 3) * (1.5 + h * (6.0 + h * (16.5 + h * (40.0 + h * (94.5 + h * (224.0 + h * (520.5 + h * (1152.0 + h * (2482.5 + h * (5632.0 + h * (13324.5 + h * (26624.0 + h * (38230.5 + h * (122880.0 + h * (726928.5 + h * (557056.0 + h * (-13230661.5 + h * (2490368.0 + h * (560069332.5 + h * (11010048.0 + h * (-25948572193.5 + h * (48234496.0 + h * (1452676184536.5 + h * (209715200.0 + h * (-95664400034365.5 + h * (905969664.0 + h * (7327814956432540.0 + h * (3892314112.0 + h * (-6.459425361709452e+17 + h * (16642998272.0 + h * (6.492408187491339e+19 + h * (70866960384.0 + h * (-7.380723366746053e+21 + h * (300647710720.0 + h * (9.422577708650279e+23 + h * (1271310319616.0 + h * (-1.3423176573208012e+26 + h * (5360119185408.0 + h * (2.1216813362745974e+28 + h * (22539988369408.0))))))))))))))))))))))))))))))))))))))));
  }
  const double PA1 = R::psigamma(0.5 * (nu + 1.0), 1.0);
  const double PB1 = R::psigamma(0.5 * nu, 1.0);
  return (1.0/4.0)*PA1 - 1.0/4.0*PB1 + (1.0/2.0)/std::pow(nu - 2, 2);
}

static inline double t2_V2(double nu) {
  if (nu >= 20.0) {
    const double h = 1.0 / nu;
    return std::pow(h, 4) * (-1.5 + h * (3.0 + h * (-21.5 + h * (21.0 + h * (-185.5 + h * (267.0 + h * (-1501.5 + h * (2721.0 + h * (-12321.5 + h * (29127.0 + h * (-104069.5 + h * (249405.0 + h * (-901129.5 + h * (2960643.0 + h * (-7935405.5 + h * (8434713.0 + h * (-70588529.5 + h * (762222975.0 + h * (-631451989.5 + h * (-24096858315.0 + h * (-5665591641.5 + h * (1469485613883.0 + h * (-50912031101.5 + h * (-95512480432239.0 + h * (-457861550785.5 + h * (7329185006684727.0 + h * (-4119232822821.5 + h * (-6.459301936838193e+17 + h * (-37066473997481.5 + h * (6.492419300919504e+19 + h * (-333569632862029.5 + h * (-7.380722366314896e+21 + h * (-3002003573362449.5 + h * (9.422577798698589e+23 + h * (-2.7017505310940404e+16 + h * (-1.342317656510326e+26 + h * (-2.431553029622236e+17 + h * (2.1216813362745974e+28 + h * (22539988369408.0)))))))))))))))))))))))))))))))))))))));
  }
  const double PA1 = R::psigamma(0.5 * (nu + 1.0), 1.0);
  const double PB1 = R::psigamma(0.5 * nu, 1.0);
  return (1.0/4.0)*PA1 - 1.0/4.0*PB1 + (-3*nu - 15)/(2*std::pow(nu, 4) - 18*std::pow(nu, 2) + 8*nu + 24) + (1.0/2.0)/std::pow(nu - 2, 2);
}

static inline double t2_V3(double nu) {
  if (nu >= 20.0) {
    const double h = 1.0 / nu;
    return std::pow(h, 4) * (-4.5 + h * (-24.0 + h * (-82.5 + h * (-240.0 + h * (-661.5 + h * (-1792.0 + h * (-4684.5 + h * (-11520.0 + h * (-27307.5 + h * (-67584.0 + h * (-173218.5 + h * (-372736.0 + h * (-573457.5 + h * (-1966080.0 + h * (-12357784.5 + h * (-10027008.0 + h * (251382568.5 + h * (-49807360.0 + h * (-11761455982.5 + h * (-242221056.0 + h * (596817160450.5 + h * (-1157627904.0 + h * (-36316904613412.5 + h * (-5452595200.0 + h * (2582938800927868.5 + h * (-25367150592.0 + h * (-2.1250663373654368e+17 + h * (-116769423360.0 + h * (2.0024218621299302e+19 + h * (-532575944704.0 + h * (-2.1424947018721418e+21 + h * (-2409476653056.0 + h * (2.5832531783611187e+23 + h * (-10823317585920.0 + h * (-3.486353752200603e+25 + h * (-48309792145408.0 + h * (5.235038863551124e+27 + h * (-214404767416320.0 + h * (-8.69889347872585e+29 + h * (-946679511515136.0))))))))))))))))))))))))))))))))))))))));
  }
  const double PA2 = R::psigamma(0.5 * (nu + 1.0), 2.0);
  const double PB2 = R::psigamma(0.5 * nu, 2.0);
  return (1.0/8.0)*PA2 - 1.0/8.0*PB2 - 1/std::pow(nu - 2, 3);
}

static inline double t2_V4(double nu) {
  if (nu >= 20.0) {
    const double h = 1.0 / nu;
    return std::pow(h, 5) * (9.0 + h * (-36.0 + h * (351.0 + h * (-1305.0 + h * (8693.0 + h * (-39402.0 + h * (213819.0 + h * (-1040499.0 + h * (5327985.0 + h * (-26505252.0 + h * (133266455.0 + h * (-665462493.0 + h * (3334101933.0 + h * (-16677900990.0 + h * (83386719411.0 + h * (-416682550503.0 + h * (2085048073385.0 + h * (-10437245145720.0 + h * (52130025018639.0 + h * (-260056733235345.0 + h * (1303287161682789.0 + h * (-6552789713368146.0 + h * (3.2582518727074092e+16 + h * (-1.6033001798896304e+17 + h * (8.14566079280641e+17 + h * (-4.285340443127665e+18 + h * (2.0364180244144456e+19 + h * (-8.179671399925876e+19 + h * (5.091047617150186e+20 + h * (-4.6880187962301177e+21 + h * (1.27276213492719e+22 + h * (1.946872085030212e+23 + h * (3.18190554516763e+23 + h * (-3.6454490317937506e+25 + h * (7.954764050109779e+24 + h * (5.195265043090143e+27 + h * (-214404767416320.0 + h * (-8.69889347872585e+29 + h * (-946679511515136.0)))))))))))))))))))))))))))))))))))))));
  }
  const double PA2 = R::psigamma(0.5 * (nu + 1.0), 2.0);
  const double PB2 = R::psigamma(0.5 * nu, 2.0);
  const double u0 = std::pow(nu, 2);
  return (1.0/8.0)*PA2 - 1.0/8.0*PB2 + (93*nu + 9*u0 + 120)/(2*std::pow(nu, 6) + 6*std::pow(nu, 5) - 38*std::pow(nu, 4) - 46*std::pow(nu, 3) - 8*nu + 228*u0 - 240) - 1/std::pow(nu - 2, 3);
}

static inline double t2_V5(double nu) {
  if (nu >= 20.0) {
    const double h = 1.0 / nu;
    return std::pow(h, 5) * (18.0 + h * (120.0 + h * (495.0 + h * (1680.0 + h * (5292.0 + h * (16128.0 + h * (46845.0 + h * (126720.0 + h * (327690.0 + h * (878592.0 + h * (2425059.0 + h * (5591040.0 + h * (9175320.0 + h * (33423360.0 + h * (222440121.0 + h * (190513152.0 + h * (-5027651370.0 + h * (1045954560.0 + h * (258752031615.0 + h * (5571084288.0 + h * (-14323611850812.0 + h * (28940697600.0 + h * (944239519948725.0 + h * (147220070400.0 + h * (-7.232228642598032e+16 + h * (735647367168.0 + h * (6.37519901209631e+18 + h * (3619852124160.0 + h * (-6.407749958815777e+20 + h * (17575006175232.0 + h * (7.284481986365282e+22 + h * (84331682856960.0 + h * (-9.299711442100027e+24 + h * (400462750679040.0 + h * (1.324814425836229e+27 + h * (1884081893670912.0 + h * (-2.09401554542045e+29 + h * (8790595464069120.0 + h * (3.6535352610648566e+31 + h * (4.070721899515085e+16))))))))))))))))))))))))))))))))))))))));
  }
  const double PA3 = R::psigamma(0.5 * (nu + 1.0), 3.0);
  const double PB3 = R::psigamma(0.5 * nu, 3.0);
  return (1.0/16.0)*PA3 - 1.0/16.0*PB3 + 3/std::pow(nu - 2, 4);
}

static inline double t2_V6(double nu) {
  if (nu >= 20.0) {
    const double h = 1.0 / nu;
    return std::pow(h, 6) * (-54.0 + h * (324.0 + h * (-4119.0 + h * (25536.0 + h * (-215472.0 + h * (1504776.0 + h * (-11089761.0 + h * (78669876.0 + h * (-560957034.0 + h * (3964217748.0 + h * (-27970619787.0 + h * (196807953912.0 + h * (-1383012798852.0 + h * (9707262418608.0 + h * (-68081417587989.0 + h * (477218778759180.0 + h * (-3343849624017150.0 + h * (2.342358259126614e+16 + h * (-1.640452232027054e+17 + h * (1.1487118100939676e+18 + h * (-8.043131181168573e+18 + h * (5.631310260678806e+19 + h * (-3.942363110349204e+20 + h * (2.7598378613088e+21 + h * (-1.9320651326450806e+22 + h * (1.3525733465893676e+23 + h * (-9.46788717136438e+23 + h * (6.627040249373051e+24 + h * (-4.639456719282255e+25 + h * (3.248388152798589e+26 + h * (-2.2733817937716025e+27 + h * (1.5904472847711494e+28 + h * (-1.1139690792783572e+29 + h * (7.811056699894509e+29 + h * (-5.458478489289088e+30 + h * (-2.09401554542045e+29 + h * (8790595464069120.0 + h * (3.6535352610648566e+31 + h * (4.070721899515085e+16)))))))))))))))))))))))))))))))))))))));
  }
  const double PA3 = R::psigamma(0.5 * (nu + 1.0), 3.0);
  const double PB3 = R::psigamma(0.5 * nu, 3.0);
  const double u0 = std::pow(nu, 2);
  const double u1 = std::pow(nu, 3);
  return (1.0/16.0)*PA3 - 1.0/16.0*PB3 + (-1239*nu - 318*u0 - 18*u1 - 1155)/(std::pow(nu, 8) + 8*std::pow(nu, 7) - 18*std::pow(nu, 6) - 160*std::pow(nu, 5) + 265*std::pow(nu, 4) - 544*nu - 1736*u0 + 888*u1 + 1680) + 3/std::pow(nu - 2, 4);
}

static inline double t2_V7(double nu) {
  if (nu >= 20.0) {
    const double h = 1.0 / nu;
    return std::pow(h, 6) * (-90.0 + h * (-720.0 + h * (-3465.0 + h * (-13440.0 + h * (-47628.0 + h * (-161280.0 + h * (-515295.0 + h * (-1520640.0 + h * (-4259970.0 + h * (-12300288.0 + h * (-36375885.0 + h * (-89456640.0 + h * (-155980440.0 + h * (-601620480.0 + h * (-4226362299.0 + h * (-3810263040.0 + h * (105580678770.0 + h * (-23011000320.0 + h * (-5951296727145.0 + h * (-133706022912.0 + h * (358090296270300.0 + h * (-752458137600.0 + h * (-2.5494467038615576e+16 + h * (-4122161971200.0 + h * (2.0973463063534292e+18 + h * (-22069421015040.0 + h * (-1.9763116937498562e+20 + h * (-115835267973120.0 + h * (2.114557486409206e+22 + h * (-597550209957888.0 + h * (-2.549568695227849e+24 + h * (-3035940582850560.0 + h * (3.44089323357701e+26 + h * (-1.521758452580352e+16 + h * (-5.166776260761293e+28 + h * (-7.536327574683648e+16 + h * (8.585463736223844e+30 + h * (-3.6920500949090304e+17 + h * (-1.5710201622578883e+33 + h * (-1.7911176357866373e+18))))))))))))))))))))))))))))))))))))))));
  }
  const double PA4 = R::psigamma(0.5 * (nu + 1.0), 4.0);
  const double PB4 = R::psigamma(0.5 * nu, 4.0);
  return (1.0/32.0)*PA4 - 1.0/32.0*PB4 - 12/std::pow(nu - 2, 5);
}

static inline double t2_V8(double nu) {
  if (nu >= 20.0) {
    const double h = 1.0 / nu;
    return std::pow(h, 5) * (6.0 + h * (-15.0 + h * (129.0 + h * (-147.0 + h * (1484.0 + h * (-2403.0 + h * (15015.0 + h * (-29931.0 + h * (147858.0 + h * (-378651.0 + h * (1456973.0 + h * (-3741075.0 + h * (14418072.0 + h * (-50330931.0 + h * (142837299.0 + h * (-160259547.0 + h * (1411770590.0 + h * (-16006682475.0 + h * (13891943769.0 + h * (554227741245.0 + h * (135974199396.0 + h * (-36737140347075.0 + h * (1323712808639.0 + h * (2578836971670453.0 + h * (12820123421994.0 + h * (-2.125463651938571e+17 + h * (123576984684645.0 + h * (2.0023836004198396e+19 + h * (1186127167919408.0 + h * (-2.1424983693034366e+21 + h * (1.1341367517309004e+16 + h * (2.5832528282102135e+23 + h * (1.0807212864104818e+17 + h * (-3.4863537855184777e+25 + h * (1.0266652018157354e+18 + h * (5.235038860390271e+27 + h * (-214404767416320.0 + h * (-8.69889347872585e+29 + h * (-946679511515136.0)))))))))))))))))))))))))))))))))))))));
  }
  const double PA2 = R::psigamma(0.5 * (nu + 1.0), 2.0);
  const double PB2 = R::psigamma(0.5 * nu, 2.0);
  const double u0 = std::pow(nu, 2);
  const double u1 = std::pow(nu, 3);
  return (1.0/8.0)*PA2 - 1.0/8.0*PB2 + (129*nu + 78*u0 + 9*u1 - 12)/(2*std::pow(nu, 7) + 4*std::pow(nu, 6) - 28*std::pow(nu, 5) - 40*std::pow(nu, 4) - 168*nu + 116*u0 + 130*u1 - 144) - 1/std::pow(nu - 2, 3);
}

static inline double t2_V9(double nu) {
  if (nu >= 20.0) {
    const double h = 1.0 / nu;
    return std::pow(h, 6) * (-30.0 + h * (90.0 + h * (-903.0 + h * (1176.0 + h * (-13356.0 + h * (24030.0 + h * (-165165.0 + h * (359172.0 + h * (-1922154.0 + h * (5301114.0 + h * (-21854595.0 + h * (59857200.0 + h * (-245107224.0 + h * (905956758.0 + h * (-2713908681.0 + h * (3205190940.0 + h * (-29647182390.0 + h * (352147014450.0 + h * (-319514706687.0 + h * (-13301465789880.0 + h * (-3399354984900.0 + h * (955165649023950.0 + h * (-35740245833253.0 + h * (-7.220743520677269e+16 + h * (-371783579237826.0 + h * (6.376390955815713e+18 + h * (-3830886525223995.0 + h * (-6.407627521343487e+20 + h * (-3.914219654134046e+16 + h * (7.284494455631684e+22 + h * (-3.969478631058151e+17 + h * (-9.299710181556769e+24 + h * (-3.998668759718783e+18 + h * (1.3248144384970216e+27 + h * (-4.003994287081368e+19 + h * (-2.09401554542045e+29 + h * (8790595464069120.0 + h * (3.6535352610648566e+31 + h * (4.070721899515085e+16)))))))))))))))))))))))))))))))))))))));
  }
  const double PA3 = R::psigamma(0.5 * (nu + 1.0), 3.0);
  const double PB3 = R::psigamma(0.5 * nu, 3.0);
  const double u0 = std::pow(nu, 2);
  const double u1 = std::pow(nu, 3);
  const double u2 = std::pow(nu, 4);
  const double u3 = std::pow(nu, 5);
  return (1.0/16.0)*PA3 - 1.0/16.0*PB3 + (-267*nu - 801*u0 - 735*u1 - 222*u2 - 18*u3 - 429)/(std::pow(nu, 10) + 4*std::pow(nu, 9) - 15*std::pow(nu, 8) - 64*std::pow(nu, 7) + 83*std::pow(nu, 6) + 864*nu - 72*u0 - 920*u1 - 173*u2 + 372*u3 + 432) + 3/std::pow(nu - 2, 4);
}


namespace d7 {

// q = z^2 ik, t = 1/(1 + q), u = q/(1 + q) and zt = z t. Every observed
// component is a sum of terms u^p t^c, u^p zt t^c and u^p z zt t^c, which
// are bounded where q^p and t^c taken apart overflow and underflow (the
// fourth derivatives were NaN from |z| of about 1e40); u is formed as
// 1/(1 + 1/q) above q = 1 and zt as 1/(1/z + z ik) above |z| = 1, where q
// itself may overflow and t underflow.
struct student_t2_TQ { double q, t, u, zt; };

inline student_t2_TQ student_t2_tq(double z, double ik) {
  student_t2_TQ T;
  T.q = z * z * ik;
  T.t = 1.0 / (1.0 + T.q);
  T.u = (T.q < 1.0) ? T.q * T.t : 1.0 / (1.0 + 1.0 / T.q);
  T.zt = (std::fabs(z) <= 1.0) ? z * T.t : 1.0 / (1.0 / z + z * ik);
  return T;
}

// D(q) = q/(1 + q) - log1p(q), with log q = 2 log|z| + log ik where q
// overflows
inline double student_t2_DQ(double z, double ik) {
  const student_t2_TQ T = student_t2_tq(z, ik);
  return std::isfinite(T.q) ? t2_D(T.q) :
    T.u - (2.0 * std::log(std::fabs(z)) + std::log(ik));
}

inline double student_t2_score_mu(double z, double s, double ik) {
  const student_t2_TQ T = student_t2_tq(z, ik);
  const double w0 = 1.0/s;
  return T.zt*w0*(3*ik + 1);
}

inline double student_t2_score_sigma(double z, double s, double ik) {
  const student_t2_TQ T = student_t2_tq(z, ik);
  const double w0 = 1.0/s;
  return w0*(-T.t + 2*T.u + z*T.zt);
}

inline double student_t2_score_nu(double z, double ik, double V0, double DQ) {
  const student_t2_TQ T = student_t2_tq(z, ik);
  return (1.0/2.0)*DQ + V0 + (3.0/2.0)*ik*T.u;
}

inline double student_t2_hess_mu_mu(double z, double s, double ik) {
  const student_t2_TQ T = student_t2_tq(z, ik);
  const double w0 = std::pow(s, -2);
  const double w1 = 3*ik;
  const double w2 = T.t*w1 - T.u;
  return T.t*w0*(3*ik*T.u - T.t - w2);
}

inline double student_t2_hess_sigma_sigma(double z, double s, double ik) {
  const student_t2_TQ T = student_t2_tq(z, ik);
  const double w0 = std::pow(s, -2);
  const double w3 = T.t*T.u;
  const double w4 = z*T.zt;
  const double w5 = std::pow(T.t, 2);
  const double w6 = std::pow(T.u, 2);
  return w0*(-3*T.t*w4 - T.u*w4 - 7*w3 + w5 - 2*w6);
}

inline double student_t2_hess_nu_nu(double z, double ik, double V1) {
  const student_t2_TQ T = student_t2_tq(z, ik);
  const double w3 = T.t*T.u;
  const double w6 = std::pow(T.u, 2);
  const double w7 = std::pow(ik, 2);
  return V1 + (1.0/2.0)*ik*w6 - 3*w3*w7 - 3.0/2.0*w6*w7;
}

inline double student_t2_expected_mu_mu(double s, double ik) {
  const double w0 = std::pow(ik, 2);
  const double w1 = 5*ik + 1;
  const double w2 = 1/(std::pow(s, 2)*w1);
  return -w2*(6*w0 + w1);
}

inline double student_t2_expected_sigma_sigma(double s, double ik) {
  const double w1 = 5*ik + 1;
  const double w2 = 1/(std::pow(s, 2)*w1);
  return -2*w2*(2*ik + 1);
}

inline double student_t2_expected_nu_nu(double V2) {
  return V2;
}

inline double student_t2_dexpected_mu_mu_mu() {
  return 0;
}

inline double student_t2_dexpected_sigma_sigma_sigma(double s, double ik) {
  const double w2 = 5*ik;
  const double w3 = w2 + 1;
  const double w4 = 1/(std::pow(s, 3)*w3);
  return 4*w4*(2*ik + 1);
}

inline double student_t2_dexpected_nu_nu_nu(double V8) {
  return V8;
}

inline void student_t2_score_curv(int k, double y, const double* th,
                                  double* out) {
  const double m = th[0], s = th[1], v = th[2];
  const double ik = 1.0 / (v - 2.0);
  const double z = (y - m) / s;
  if (k == 0) {
    out[0] = student_t2_score_mu(z, s, ik);
    out[1] = student_t2_hess_mu_mu(z, s, ik);
  } else if (k == 1) {
    out[0] = student_t2_score_sigma(z, s, ik);
    out[1] = student_t2_hess_sigma_sigma(z, s, ik);
  } else {
    out[0] = student_t2_score_nu(z, ik, t2_V0(v), student_t2_DQ(z, ik));
    out[1] = student_t2_hess_nu_nu(z, ik, t2_V1(v));
  }
}

inline void student_t2_info_dinfo(int k, double y, const double* th,
                                  double* out) {
  const double s = th[1], v = th[2];
  const double ik = 1.0 / (v - 2.0);
  if (k == 0) {
    out[0] = student_t2_expected_mu_mu(s, ik);
    out[1] = student_t2_dexpected_mu_mu_mu();
  } else if (k == 1) {
    out[0] = student_t2_expected_sigma_sigma(s, ik);
    out[1] = student_t2_dexpected_sigma_sigma_sigma(s, ik);
  } else {
    out[0] = student_t2_expected_nu_nu(t2_V2(v));
    out[1] = student_t2_dexpected_nu_nu_nu(t2_V8(v));
  }
}

} // namespace d7

#endif
