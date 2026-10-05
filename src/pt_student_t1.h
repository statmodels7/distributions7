#ifndef D7_PT_STUDENT_T1_H
#define D7_PT_STUDENT_T1_H

#include <Rcpp.h>
#include <cmath>

// The Student t in its location mu, scale s and degrees of freedom nu:
// one function per component and order for the quantities the scalar
// registry reads, lifted from the generated kernels in student_t1.cpp (which
// records the derivation) and called by them and by the registry
// (d7_ccallable.cpp). With ik = 1/nu and z = (y - mu)/s, the caller
// passes z, ik, D(q) for q = z^2 ik, and the quantities of nu alone
// t1_V*(nu), each of which carries a polygamma.

// D(q) = q/(1+q) - log1p(q). The direct form cancels as q -> 0, where D is
// -q^2/2. With w = q/(2+q), q/(1+q) = 2w/(1+w) and log1p(q) = 2 atanh(w), so
//   D = -2w^2/(1+w) - 2 sum_{j>=1} w^(2j+1)/(2j+1),
// a sum of terms of one sign; below q = 0.5 (w <= 0.2) twelve terms reach
// the last bit, and above it the direct form loses at most a few.
static inline double t1_D(double q) {
  if (q < 0.5) {
    const double w = q / (2.0 + q), w2 = w * w;
    double acc = 1.0 / 25.0;
    for (int j = 11; j >= 1; j--) acc = 1.0 / (2.0 * j + 1.0) + w2 * acc;
    return -2.0 * w2 / (1.0 + w) - 2.0 * w * w2 * acc;
  }
  return q / (1.0 + q) - std::log1p(q);
}

static inline double t1_V0(double nu) {
  if (nu >= 20.0) {
    const double h = 1.0 / nu;
    return std::pow(h, 2) * (0.25 + h * (0.0 + h * (-0.125 + h * (0.0 + h * (0.25 + h * (0.0 + h * (-1.0625 + h * (0.0 + h * (7.75 + h * (0.0 + h * (-86.375 + h * (0.0 + h * (1365.25 + h * (0.0 + h * (-29049.03125 + h * (0.0 + h * (800572.75 + h * (0.0 + h * (-27741322.625 + h * (0.0 + h * (1180529130.25 + h * (0.0 + h * (-60523980051.6875 + h * (0.0 + h * (3679416778537.75 + h * (0.0 + h * (-261707609906583.88 + h * (0.0 + h * (2.1531418140800296e+16 + h * (0.0 + h * (-2.0288775575173015e+18 + h * (0.0 + h * (2.1708009902623772e+20 + h * (0.0 + h * (-2.6173826968455817e+22 + h * (0.0 + h * (3.532414887686388e+24 + h * (0.0 + h * (-5.3042033406864905e+26)))))))))))))))))))))))))))))))))))))));
  }
  const double PA0 = R::psigamma(0.5 * (nu + 1.0), 0.0);
  const double PB0 = R::psigamma(0.5 * nu, 0.0);
  return (1.0/2.0)*PA0 - 1.0/2.0*PB0 - (1.0/2.0)/nu;
}

static inline double t1_V1(double nu) {
  if (nu >= 20.0) {
    const double h = 1.0 / nu;
    return std::pow(h, 3) * (-0.5 + h * (0.0 + h * (0.5 + h * (0.0 + h * (-1.5 + h * (0.0 + h * (8.5 + h * (0.0 + h * (-77.5 + h * (0.0 + h * (1036.5 + h * (0.0 + h * (-19113.5 + h * (0.0 + h * (464784.5 + h * (0.0 + h * (-14410309.5 + h * (0.0 + h * (554826452.5 + h * (0.0 + h * (-25971640865.5 + h * (0.0 + h * (1452575521240.5 + h * (0.0 + h * (-95664836241981.5 + h * (0.0 + h * (7327813077384348.0 + h * (0.0 + h * (-6.459425442240088e+17 + h * (0.0 + h * (6.492408184055365e+19 + h * (0.0 + h * (-7.380723366892082e+21 + h * (0.0 + h * (9.422577708644094e+23 + h * (0.0 + h * (-1.3423176573208274e+26 + h * (0.0 + h * (2.1216813362745965e+28)))))))))))))))))))))))))))))))))))))));
  }
  const double PA1 = R::psigamma(0.5 * (nu + 1.0), 1.0);
  const double PB1 = R::psigamma(0.5 * nu, 1.0);
  return (1.0/4.0)*PA1 - 1.0/4.0*PB1 + (1.0/2.0)/std::pow(nu, 2);
}

static inline double t1_V2(double nu) {
  if (nu >= 20.0) {
    const double h = 1.0 / nu;
    return std::pow(h, 4) * (-3.5 + h * (13.0 + h * (-39.5 + h * (119.0 + h * (-363.5 + h * (1101.0 + h * (-3279.5 + h * (9763.0 + h * (-29523.5 + h * (89609.0 + h * (-265719.5 + h * (778047.0 + h * (-2391483.5 + h * (7639237.0 + h * (-21523359.5 + h * (50159771.0 + h * (-193710243.5 + h * (1135957185.0 + h * (-1743392199.5 + h * (-20741464265.0 + h * (-15690529803.5 + h * (1499647110653.0 + h * (-141214768239.5 + h * (-95241191937261.0 + h * (-1270932914163.5 + h * (7331625876126841.0 + h * (-11438396227479.5 + h * (-6.459082290353265e+17 + h * (-102945566047323.5 + h * (6.492439067725179e+19 + h * (-926510094425919.5 + h * (-7.380720587361798e+21 + h * (-8338590849833284.0 + h * (9.422577958801819e+23 + h * (-7.504731764849955e+16 + h * (-1.3423176550694079e+26 + h * (-6.75425858836496e+17 + h * (2.1216813362745965e+28))))))))))))))))))))))))))))))))))))));
  }
  const double PA1 = R::psigamma(0.5 * (nu + 1.0), 1.0);
  const double PB1 = R::psigamma(0.5 * nu, 1.0);
  const double u0 = std::pow(nu, 2);
  return (1.0/4.0)*PA1 - 1.0/4.0*PB1 + (nu - 3)/(2*std::pow(nu, 4) + 8*std::pow(nu, 3) + 6*u0) + (1.0/2.0)/u0;
}

static inline double t1_V3(double nu) {
  if (nu >= 20.0) {
    const double h = 1.0 / nu;
    return std::pow(h, 4) * (1.5 + h * (0.0 + h * (-2.5 + h * (0.0 + h * (10.5 + h * (0.0 + h * (-76.5 + h * (0.0 + h * (852.5 + h * (0.0 + h * (-13474.5 + h * (0.0 + h * (286702.5 + h * (0.0 + h * (-7901336.5 + h * (0.0 + h * (273795880.5 + h * (0.0 + h * (-11651355502.5 + h * (0.0 + h * (597347739906.5 + h * (0.0 + h * (-36314388031012.5 + h * (0.0 + h * (2582950578533500.5 + h * (0.0 + h * (-2.125065792441461e+17 + h * (0.0 + h * (2.0024218870944276e+19 + h * (0.0 + h * (-2.1424947007382705e+21 + h * (0.0 + h * (2.5832531784122288e+23 + h * (0.0 + h * (-3.4863537521983147e+25 + h * (0.0 + h * (5.235038863551227e+27 + h * (0.0 + h * (-8.698893478725845e+29)))))))))))))))))))))))))))))))))))))));
  }
  const double PA2 = R::psigamma(0.5 * (nu + 1.0), 2.0);
  const double PB2 = R::psigamma(0.5 * nu, 2.0);
  return (1.0/8.0)*PA2 - 1.0/8.0*PB2 - 1/std::pow(nu, 3);
}

static inline double t1_V4(double nu) {
  if (nu >= 20.0) {
    const double h = 1.0 / nu;
    return std::pow(h, 5) * (21.0 + h * (-142.0 + h * (795.0 + h * (-4251.0 + h * (22161.0 + h * (-113436.0 + h * (574455.0 + h * (-2894389.0 + h * (14545101.0 + h * (-72945654.0 + h * (365280915.0 + h * (-1827977919.0 + h * (9146903241.0 + h * (-45759157936.0 + h * (228806504175.0 + h * (-1143909388521.0 + h * (5721367912581.0 + h * (-28619846890122.0 + h * (143045045588235.0 + h * (-714640083946675.0 + h * (3576223769669121.0 + h * (-1.7917543070085252e+16 + h * (8.940647291139709e+16 + h * (-4.444504024818296e+17 + h * (2.2351697308119488e+18 + h * (-1.138836412983429e+19 + h * (5.587931444254192e+19 + h * (-2.593724334105389e+20 + h * (1.3969835016137367e+21 + h * (-9.127412929425916e+21 + h * (3.4924593305295113e+22 + h * (8.370234482917664e+22 + h * (8.731148845169431e+23 + h * (-3.922911200293799e+25 + h * (2.1827872579884668e+25 + h * (5.125899500126472e+27 + h * (0.0 + h * (-8.698893478725845e+29))))))))))))))))))))))))))))))))))))));
  }
  const double PA2 = R::psigamma(0.5 * (nu + 1.0), 2.0);
  const double PB2 = R::psigamma(0.5 * nu, 2.0);
  const double u0 = std::pow(nu, 3);
  return (1.0/8.0)*PA2 - 1.0/8.0*PB2 + (-3*std::pow(nu, 2) + 15*nu + 30)/(2*std::pow(nu, 6) + 18*std::pow(nu, 5) + 46*std::pow(nu, 4) + 30*u0) - 1/u0;
}

static inline double t1_V5(double nu) {
  if (nu >= 20.0) {
    const double h = 1.0 / nu;
    return std::pow(h, 5) * (-6.0 + h * (0.0 + h * (15.0 + h * (0.0 + h * (-84.0 + h * (0.0 + h * (765.0 + h * (0.0 + h * (-10230.0 + h * (0.0 + h * (188643.0 + h * (0.0 + h * (-4587240.0 + h * (0.0 + h * (142224057.0 + h * (0.0 + h * (-5475917610.0 + h * (0.0 + h * (256329821055.0 + h * (0.0 + h * (-14336345757756.0 + h * (0.0 + h * (944174088806325.0 + h * (0.0 + h * (-7.232261619893802e+16 + h * (0.0 + h * (6.375197377324383e+18 + h * (0.0 + h * (-6.407750038702168e+20 + h * (0.0 + h * (7.28448198251012e+22 + h * (0.0 + h * (-9.299711442284024e+24 + h * (0.0 + h * (1.3248144258353596e+27 + h * (0.0 + h * (-2.0940155454204906e+29 + h * (0.0 + h * (3.653535261064855e+31)))))))))))))))))))))))))))))))))))))));
  }
  const double PA3 = R::psigamma(0.5 * (nu + 1.0), 3.0);
  const double PB3 = R::psigamma(0.5 * nu, 3.0);
  return (1.0/16.0)*PA3 - 1.0/16.0*PB3 + 3/std::pow(nu, 4);
}

static inline double t1_V6(double nu) {
  if (nu >= 20.0) {
    const double h = 1.0 / nu;
    return std::pow(h, 6) * (-126.0 + h * (1236.0 + h * (-10071.0 + h * (77592.0 + h * (-578376.0 + h * (4218936.0 + h * (-30363921.0 + h * (216687996.0 + h * (-1537543026.0 + h * (10865971764.0 + h * (-76574505771.0 + h * (538584415536.0 + h * (-3782952359676.0 + h * (26544931528128.0 + h * (-186134067315621.0 + h * (1304535490993716.0 + h * (-9139798601188326.0 + h * (6.401890376508208e+16 + h * (-4.483308123702335e+17 + h * (3.1393027344289864e+18 + h * (-2.198022636403699e+19 + h * (1.5388756290851134e+20 + h * (-1.0773315016010599e+21 + h * (7.541874039479583e+21 + h * (-5.279675378487585e+22 + h * (3.695992979298258e+23 + h * (-2.587228690280599e+24 + h * (1.8110351212613238e+25 + h * (-1.2677889969124117e+26 + h * (8.87534921547091e+26 + h * (-6.212283431530661e+27 + h * (4.347692878146207e+28 + h * (-3.044048218112514e+29 + h * (2.1321646789091162e+30 + h * (-1.491590961040532e+31 + h * (-2.0940155454204906e+29 + h * (0.0 + h * (3.653535261064855e+31))))))))))))))))))))))))))))))))))))));
  }
  const double PA3 = R::psigamma(0.5 * (nu + 1.0), 3.0);
  const double PB3 = R::psigamma(0.5 * nu, 3.0);
  const double u0 = std::pow(nu, 4);
  return (1.0/16.0)*PA3 - 1.0/16.0*PB3 + (6*std::pow(nu, 3) - 30*std::pow(nu, 2) - 279*nu - 315)/(std::pow(nu, 8) + 16*std::pow(nu, 7) + 86*std::pow(nu, 6) + 176*std::pow(nu, 5) + 105*u0) + 3/u0;
}

static inline double t1_V7(double nu) {
  if (nu >= 20.0) {
    const double h = 1.0 / nu;
    return std::pow(h, 6) * (30.0 + h * (0.0 + h * (-105.0 + h * (0.0 + h * (756.0 + h * (0.0 + h * (-8415.0 + h * (0.0 + h * (132990.0 + h * (0.0 + h * (-2829645.0 + h * (0.0 + h * (77983080.0 + h * (0.0 + h * (-2702257083.0 + h * (0.0 + h * (114994269810.0 + h * (0.0 + h * (-5895585884265.0 + h * (0.0 + h * (358408643943900.0 + h * (0.0 + h * (-2.5492700397770776e+16 + h * (0.0 + h * (2.0973558697692024e+18 + h * (0.0 + h * (-1.9763111869705588e+20 + h * (0.0 + h * (2.1145575127717155e+22 + h * (0.0 + h * (-2.549568693878542e+24 + h * (0.0 + h * (3.440893233645089e+26 + h * (0.0 + h * (-5.166776260757902e+28 + h * (0.0 + h * (8.585463736224012e+30 + h * (0.0 + h * (-1.5710201622578874e+33)))))))))))))))))))))))))))))))))))))));
  }
  const double PA4 = R::psigamma(0.5 * (nu + 1.0), 4.0);
  const double PB4 = R::psigamma(0.5 * nu, 4.0);
  return (1.0/32.0)*PA4 - 1.0/32.0*PB4 - 12/std::pow(nu, 5);
}

static inline double t1_V8(double nu) {
  if (nu >= 20.0) {
    const double h = 1.0 / nu;
    return std::pow(h, 5) * (14.0 + h * (-65.0 + h * (237.0 + h * (-833.0 + h * (2908.0 + h * (-9909.0 + h * (32795.0 + h * (-107393.0 + h * (354282.0 + h * (-1164917.0 + h * (3720073.0 + h * (-11670705.0 + h * (38263736.0 + h * (-129867029.0 + h * (387420471.0 + h * (-953035649.0 + h * (3874204870.0 + h * (-23855100885.0 + h * (38354628389.0 + h * (477053678095.0 + h * (376572715284.0 + h * (-37491177766325.0 + h * (3671583974227.0 + h * (2571512182306047.0 + h * (35586121596578.0 + h * (-2.126171504076784e+17 + h * (343151886824385.0 + h * (2.002315510009512e+19 + h * (3294258113514352.0 + h * (-2.1425048923493093e+21 + h * (3.1501343210481264e+16 + h * (2.5832522055766295e+23 + h * (3.001892705939982e+17 + h * (-3.486353844756673e+25 + h * (2.8517980706429834e+18 + h * (5.235038854770691e+27 + h * (0.0 + h * (-8.698893478725845e+29))))))))))))))))))))))))))))))))))))));
  }
  const double PA2 = R::psigamma(0.5 * (nu + 1.0), 2.0);
  const double PB2 = R::psigamma(0.5 * nu, 2.0);
  const double u0 = std::pow(nu, 3);
  return (1.0/8.0)*PA2 - 1.0/8.0*PB2 + (4*std::pow(nu, 2) + 33*nu - 3*u0 + 18)/(2*std::pow(nu, 7) + 16*std::pow(nu, 6) + 44*std::pow(nu, 5) + 48*std::pow(nu, 4) + 18*u0) - 1/u0;
}

static inline double t1_V9(double nu) {
  if (nu >= 20.0) {
    const double h = 1.0 / nu;
    return std::pow(h, 6) * (-70.0 + h * (390.0 + h * (-1659.0 + h * (6664.0 + h * (-26172.0 + h * (99090.0 + h * (-360745.0 + h * (1288716.0 + h * (-4605666.0 + h * (16308838.0 + h * (-55801095.0 + h * (186731280.0 + h * (-650483512.0 + h * (2337606522.0 + h * (-7360988949.0 + h * (19060712980.0 + h * (-81358302270.0 + h * (524812219470.0 + h * (-882156452947.0 + h * (-11449288274280.0 + h * (-9414317882100.0 + h * (974770621924450.0 + h * (-99132767304129.0 + h * (-7.200234110456931e+16 + h * (-1031997526300762.0 + h * (6.378514512230352e+18 + h * (-1.0637708491555936e+16 + h * (-6.407409632030438e+20 + h * (-1.0871051774597362e+17 + h * (7.284516633987651e+22 + h * (-1.1025470123668442e+18 + h * (-9.299707940075866e+24 + h * (-1.1107003011977933e+19 + h * (1.324814461007536e+27 + h * (-1.1122012475507635e+20 + h * (-2.0940155454204906e+29 + h * (0.0 + h * (3.653535261064855e+31))))))))))))))))))))))))))))))))))))));
  }
  const double PA3 = R::psigamma(0.5 * (nu + 1.0), 3.0);
  const double PB3 = R::psigamma(0.5 * nu, 3.0);
  const double u0 = std::pow(nu, 4);
  const double u1 = std::pow(nu, 5);
  return (1.0/16.0)*PA3 - 1.0/16.0*PB3 + (-123*std::pow(nu, 3) - 333*std::pow(nu, 2) - 279*nu + 2*u0 + 6*u1 - 81)/(std::pow(nu, 10) + 12*std::pow(nu, 9) + 57*std::pow(nu, 8) + 136*std::pow(nu, 7) + 171*std::pow(nu, 6) + 27*u0 + 108*u1) + 3/u0;
}


namespace d7 {

// q = z^2 ik, t = 1/(1 + q), u = q/(1 + q) and zt = z t. Every observed
// component is a sum of terms u^p t^c, u^p zt t^c and u^p z zt t^c, which
// are bounded where q^p and t^c taken apart overflow and underflow (the
// fourth derivatives were NaN from |z| of about 1e40); u is formed as
// 1/(1 + 1/q) above q = 1 and zt as 1/(1/z + z ik) above |z| = 1, where q
// itself may overflow and t underflow.
struct student_t1_TQ { double q, t, u, zt; };

inline student_t1_TQ student_t1_tq(double z, double ik) {
  student_t1_TQ T;
  T.q = z * z * ik;
  T.t = 1.0 / (1.0 + T.q);
  T.u = (T.q < 1.0) ? T.q * T.t : 1.0 / (1.0 + 1.0 / T.q);
  T.zt = (std::fabs(z) <= 1.0) ? z * T.t : 1.0 / (1.0 / z + z * ik);
  return T;
}

// D(q) = q/(1 + q) - log1p(q), with log q = 2 log|z| + log ik where q
// overflows
inline double student_t1_DQ(double z, double ik) {
  const student_t1_TQ T = student_t1_tq(z, ik);
  return std::isfinite(T.q) ? t1_D(T.q) :
    T.u - (2.0 * std::log(std::fabs(z)) + std::log(ik));
}

inline double student_t1_score_mu(double z, double s, double ik) {
  const student_t1_TQ T = student_t1_tq(z, ik);
  const double w0 = 1.0/s;
  return T.zt*w0*(ik + 1);
}

inline double student_t1_score_sigma(double z, double s, double ik) {
  const student_t1_TQ T = student_t1_tq(z, ik);
  const double w0 = 1.0/s;
  return w0*(-T.t + z*T.zt);
}

inline double student_t1_score_nu(double z, double ik, double V0, double DQ) {
  const student_t1_TQ T = student_t1_tq(z, ik);
  return (1.0/2.0)*DQ + V0 + (1.0/2.0)*ik*T.u;
}

inline double student_t1_hess_mu_mu(double z, double s, double ik) {
  const student_t1_TQ T = student_t1_tq(z, ik);
  const double w0 = std::pow(s, -2);
  const double w1 = ik*T.t;
  const double w2 = -T.u + w1;
  return T.t*w0*(ik*T.u - T.t - w2);
}

inline double student_t1_hess_sigma_sigma(double z, double s, double ik) {
  const student_t1_TQ T = student_t1_tq(z, ik);
  const double w0 = std::pow(s, -2);
  const double w3 = T.t*T.u;
  const double w4 = z*T.zt;
  const double w5 = std::pow(T.t, 2);
  return w0*(-3*T.t*w4 - T.u*w4 - w3 + w5);
}

inline double student_t1_hess_nu_nu(double z, double ik, double V1) {
  const student_t1_TQ T = student_t1_tq(z, ik);
  const double w3 = T.t*T.u;
  const double w6 = (1.0/2.0)*std::pow(T.u, 2);
  const double w7 = std::pow(ik, 2);
  return V1 + ik*w6 - w3*w7 - w6*w7;
}

inline double student_t1_expected_mu_mu(double s, double ik) {
  const double w0 = 1/(std::pow(s, 2)*(3*ik + 1));
  return -w0*(ik + 1);
}

inline double student_t1_expected_sigma_sigma(double s, double ik) {
  const double w0 = 1/(std::pow(s, 2)*(3*ik + 1));
  return -2*w0;
}

inline double student_t1_expected_nu_nu(double V2) {
  return V2;
}

inline double student_t1_dexpected_mu_mu_mu() {
  return 0;
}

inline double student_t1_dexpected_sigma_sigma_sigma(double s, double ik) {
  const double w0 = 1/(std::pow(s, 3)*(3*ik + 1));
  return 4*w0;
}

inline double student_t1_dexpected_nu_nu_nu(double V8) {
  return V8;
}

inline void student_t1_score_curv(int k, double y, const double* th,
                                  double* out) {
  const double m = th[0], s = th[1], v = th[2];
  const double ik = 1.0 / (v - 0.0);
  const double z = (y - m) / s;
  if (k == 0) {
    out[0] = student_t1_score_mu(z, s, ik);
    out[1] = student_t1_hess_mu_mu(z, s, ik);
  } else if (k == 1) {
    out[0] = student_t1_score_sigma(z, s, ik);
    out[1] = student_t1_hess_sigma_sigma(z, s, ik);
  } else {
    out[0] = student_t1_score_nu(z, ik, t1_V0(v), student_t1_DQ(z, ik));
    out[1] = student_t1_hess_nu_nu(z, ik, t1_V1(v));
  }
}

inline void student_t1_info_dinfo(int k, double y, const double* th,
                                  double* out) {
  const double s = th[1], v = th[2];
  const double ik = 1.0 / (v - 0.0);
  if (k == 0) {
    out[0] = student_t1_expected_mu_mu(s, ik);
    out[1] = student_t1_dexpected_mu_mu_mu();
  } else if (k == 1) {
    out[0] = student_t1_expected_sigma_sigma(s, ik);
    out[1] = student_t1_dexpected_sigma_sigma_sigma(s, ik);
  } else {
    out[0] = student_t1_expected_nu_nu(t1_V2(v));
    out[1] = student_t1_dexpected_nu_nu_nu(t1_V8(v));
  }
}

} // namespace d7

#endif
