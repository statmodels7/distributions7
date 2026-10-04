#include <Rcpp.h>
#include <cmath>
#include "d7_par.h"
using namespace Rcpp;

// The Student t in its location mu, scale s and degrees of freedom nu > 0.
// With z = (y - mu)/s, k = nu, q = z^2/k and t = 1/(1 + q), the
// log-density is
//   c(nu) - log s - log(pi)/2 - (k + 1)/2 log(1 + q),
//   c(nu) = lgamma((nu+1)/2) - lgamma(nu/2) - log(nu)/2.
// Every component below is a closed form derived offline with sympy
// (stabilita/gen_student_t2.py, T_FAMILY=t1), one exported function per
// order and per surface, and is checked numerically in the tests.
//
// Every derivative in nu vanishes as nu grows and is a difference of terms
// agreeing to leading order. The data part is cancelled symbolically and
// written in q, z, 1/k and t, which neither cancel nor overflow at the nu
// the link can produce; the score's logarithm enters as
// D(q) = q/(1+q) - log1p(q). Each quantity of nu alone that carries a
// polygamma, t1_V*, has two branches: the direct form below nu = 20, and
// above it the asymptotic series in h = 1/nu, exact to h^40, from
// Stirling's series. The expectations are closed, t being Beta(nu/2, 1/2)
// under the model.

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

// order 1 of log f in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t1_gradient_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu(n);
  NumericVector o_sigma(n);
  NumericVector o_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, V0; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double V0 = t1_V0(v);
    return Prm{m, s, ik, V0};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double V0 = P.V0; (void) V0;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    (void) t; (void) q;
    const double DQ = t1_D(q);
    const double w0 = t/s;
    o_mu[i] = w0*z*(ik + 1);
    o_sigma[i] = w0*(std::pow(z, 2) - 1);
    o_nu[i] = (1.0/2.0)*DQ + V0 + (1.0/2.0)*ik*q*t;
  });
  return List::create(Named("mu") = o_mu, Named("sigma") = o_sigma, Named("nu") = o_nu);
}

// order 2 of log f in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t1_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu(n);
  NumericVector o_sigma_sigma(n);
  NumericVector o_nu_nu(n);
  NumericVector o_mu_sigma(n);
  NumericVector o_mu_nu(n);
  NumericVector o_sigma_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, V1; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double V1 = t1_V1(v);
    return Prm{m, s, ik, V1};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double V1 = P.V1; (void) V1;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    (void) t; (void) q;
    const double w0 = q - 1;
    const double w1 = std::pow(t, 2);
    const double w2 = w1/std::pow(s, 2);
    const double w3 = std::pow(z, 2);
    const double w4 = (1.0/2.0)*std::pow(q, 2);
    const double w5 = ik*w1;
    const double w6 = std::pow(ik, 2)*w1;
    const double w7 = (-ik + q)/s;
    o_mu_mu[i] = w2*(ik*q - ik + w0);
    o_sigma_sigma[i] = w2*(-q*w3 - w0 - 3*w3);
    o_nu_nu[i] = V1 - q*w6 + w4*w5 - w4*w6;
    o_mu_sigma[i] = -2*w2*z*(ik + 1);
    o_mu_nu[i] = w5*w7*z;
    o_sigma_nu[i] = q*w1*w7;
  });
  return List::create(Named("mu_mu") = o_mu_mu, Named("sigma_sigma") = o_sigma_sigma, Named("nu_nu") = o_nu_nu, Named("mu_sigma") = o_mu_sigma, Named("mu_nu") = o_mu_nu, Named("sigma_nu") = o_sigma_nu);
}

// order 2, expected
// [[Rcpp::export]]
List student_t1_expected_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu(n);
  NumericVector o_sigma_sigma(n);
  NumericVector o_nu_nu(n);
  NumericVector o_mu_sigma(n);
  NumericVector o_mu_nu(n);
  NumericVector o_sigma_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, V2, w0, w1; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double V2 = t1_V2(v);
    const double w0 = 1/(std::pow(s, 2)*(3*ik + 1));
    const double w1 = std::pow(ik, 2);
    return Prm{m, s, ik, V2, w0, w1};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double V2 = P.V2; (void) V2;
    const double w0 = P.w0; (void) w0;
    const double w1 = P.w1; (void) w1;
    o_mu_mu[i] = -w0*(ik + 1);
    o_sigma_sigma[i] = -2*w0;
    o_nu_nu[i] = V2;
    o_mu_sigma[i] = 0;
    o_mu_nu[i] = 0;
    o_sigma_nu[i] = 2*w1/(s*(4*ik + 3*w1 + 1));
  });
  return List::create(Named("mu_mu") = o_mu_mu, Named("sigma_sigma") = o_sigma_sigma, Named("nu_nu") = o_nu_nu, Named("mu_sigma") = o_mu_sigma, Named("mu_nu") = o_mu_nu, Named("sigma_nu") = o_sigma_nu);
}

// order 3 of log f in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t1_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu_mu(n);
  NumericVector o_mu_mu_sigma(n);
  NumericVector o_mu_mu_nu(n);
  NumericVector o_mu_sigma_sigma(n);
  NumericVector o_mu_sigma_nu(n);
  NumericVector o_mu_nu_nu(n);
  NumericVector o_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_nu(n);
  NumericVector o_sigma_nu_nu(n);
  NumericVector o_nu_nu_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, V3, w1, w10, w13, w16; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double V3 = t1_V3(v);
    const double w1 = 3*ik;
    const double w10 = -ik;
    const double w13 = std::pow(s, -2);
    const double w16 = std::pow(ik, 2);
    return Prm{m, s, ik, V3, w1, w10, w13, w16};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double V3 = P.V3; (void) V3;
    const double w1 = P.w1; (void) w1;
    const double w10 = P.w10; (void) w10;
    const double w13 = P.w13; (void) w13;
    const double w16 = P.w16; (void) w16;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    (void) t; (void) q;
    const double w0 = ik*q;
    const double w2 = w0 - w1;
    const double w3 = q + w2 - 3;
    const double w4 = std::pow(t, 3);
    const double w5 = 2*w4;
    const double w6 = w5/std::pow(s, 3);
    const double w7 = w6*z;
    const double w8 = 3*q;
    const double w9 = w8 - 1;
    const double w11 = q*w1 + w10;
    const double w12 = std::pow(q, 2);
    const double w14 = w13*w4;
    const double w15 = w5*z;
    const double w17 = (ik - q)/s;
    const double w18 = std::pow(z, 2);
    const double w19 = std::pow(ik, 3)*w4;
    const double w20 = std::pow(q, 3);
    const double w21 = w16*w4;
    o_mu_mu_mu[i] = ik*w3*w7;
    o_mu_mu_sigma[i] = w6*(-w11 - w9);
    o_mu_mu_nu[i] = ik*w14*(-w11 + w12 - w8);
    o_mu_sigma_sigma[i] = -w3*w7;
    o_mu_sigma_nu[i] = ik*w13*w15*(-2*q - w0 - w10);
    o_mu_nu_nu[i] = w15*w16*w17;
    o_sigma_sigma_sigma[i] = w6*(w12*w18 + w18*w8 + 6*w18 + w9);
    o_sigma_sigma_nu[i] = q*w14*(-5*q - w12 - w2);
    o_sigma_nu_nu[i] = w0*w17*w5;
    o_nu_nu_nu[i] = V3 + 3*w12*w19 - 3.0/2.0*w12*w21 + w19*w20 + w19*w8 - 1.0/2.0*w20*w21;
  });
  return List::create(Named("mu_mu_mu") = o_mu_mu_mu, Named("mu_mu_sigma") = o_mu_mu_sigma, Named("mu_mu_nu") = o_mu_mu_nu, Named("mu_sigma_sigma") = o_mu_sigma_sigma, Named("mu_sigma_nu") = o_mu_sigma_nu, Named("mu_nu_nu") = o_mu_nu_nu, Named("sigma_sigma_sigma") = o_sigma_sigma_sigma, Named("sigma_sigma_nu") = o_sigma_sigma_nu, Named("sigma_nu_nu") = o_sigma_nu_nu, Named("nu_nu_nu") = o_nu_nu_nu);
}

// order 3, expected
// [[Rcpp::export]]
List student_t1_deriv3_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu_mu(n);
  NumericVector o_mu_mu_sigma(n);
  NumericVector o_mu_mu_nu(n);
  NumericVector o_mu_sigma_sigma(n);
  NumericVector o_mu_sigma_nu(n);
  NumericVector o_mu_nu_nu(n);
  NumericVector o_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_nu(n);
  NumericVector o_sigma_nu_nu(n);
  NumericVector o_nu_nu_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, V4, w0, w1, w2, w3, w4, w5, w6; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double V4 = t1_V4(v);
    const double w0 = std::pow(ik, 2);
    const double w1 = 2*w0;
    const double w2 = 1.0/(8*ik + 15*w0 + 1);
    const double w3 = 2*w2/std::pow(s, 3);
    const double w4 = ik - 1;
    const double w5 = w2/std::pow(s, 2);
    const double w6 = std::pow(ik, 3);
    return Prm{m, s, ik, V4, w0, w1, w2, w3, w4, w5, w6};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double V4 = P.V4; (void) V4;
    const double w0 = P.w0; (void) w0;
    const double w1 = P.w1; (void) w1;
    const double w2 = P.w2; (void) w2;
    const double w3 = P.w3; (void) w3;
    const double w4 = P.w4; (void) w4;
    const double w5 = P.w5; (void) w5;
    const double w6 = P.w6; (void) w6;
    o_mu_mu_mu[i] = 0;
    o_mu_mu_sigma[i] = w3*(3*ik + w1 + 1);
    o_mu_mu_nu[i] = w1*w4*w5;
    o_mu_sigma_sigma[i] = 0;
    o_mu_sigma_nu[i] = 0;
    o_mu_nu_nu[i] = 0;
    o_sigma_sigma_sigma[i] = w3*(13*ik + 5);
    o_sigma_sigma_nu[i] = -12*w0*w5;
    o_sigma_nu_nu[i] = 4*w4*w6/(s*(9*ik + 23*w0 + 15*w6 + 1));
    o_nu_nu_nu[i] = V4;
  });
  return List::create(Named("mu_mu_mu") = o_mu_mu_mu, Named("mu_mu_sigma") = o_mu_mu_sigma, Named("mu_mu_nu") = o_mu_mu_nu, Named("mu_sigma_sigma") = o_mu_sigma_sigma, Named("mu_sigma_nu") = o_mu_sigma_nu, Named("mu_nu_nu") = o_mu_nu_nu, Named("sigma_sigma_sigma") = o_sigma_sigma_sigma, Named("sigma_sigma_nu") = o_sigma_sigma_nu, Named("sigma_nu_nu") = o_sigma_nu_nu, Named("nu_nu_nu") = o_nu_nu_nu);
}

// order 4 of log f in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t1_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu_mu_mu(n);
  NumericVector o_mu_mu_mu_sigma(n);
  NumericVector o_mu_mu_mu_nu(n);
  NumericVector o_mu_mu_sigma_sigma(n);
  NumericVector o_mu_mu_sigma_nu(n);
  NumericVector o_mu_mu_nu_nu(n);
  NumericVector o_mu_sigma_sigma_sigma(n);
  NumericVector o_mu_sigma_sigma_nu(n);
  NumericVector o_mu_sigma_nu_nu(n);
  NumericVector o_mu_nu_nu_nu(n);
  NumericVector o_sigma_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_sigma_nu(n);
  NumericVector o_sigma_sigma_nu_nu(n);
  NumericVector o_sigma_nu_nu_nu(n);
  NumericVector o_nu_nu_nu_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, V5, w7, w10, w14, w15, w22, w23, w24, w28, w35; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double V5 = t1_V5(v);
    const double w7 = 6*ik;
    const double w10 = -ik;
    const double w14 = std::pow(ik, 2);
    const double w15 = std::pow(s, -3);
    const double w22 = std::pow(s, -2);
    const double w23 = w14*w22;
    const double w24 = 3*ik;
    const double w28 = std::pow(ik, 3);
    const double w35 = std::pow(ik, 4);
    return Prm{m, s, ik, V5, w7, w10, w14, w15, w22, w23, w24, w28, w35};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double V5 = P.V5; (void) V5;
    const double w7 = P.w7; (void) w7;
    const double w10 = P.w10; (void) w10;
    const double w14 = P.w14; (void) w14;
    const double w15 = P.w15; (void) w15;
    const double w22 = P.w22; (void) w22;
    const double w23 = P.w23; (void) w23;
    const double w24 = P.w24; (void) w24;
    const double w28 = P.w28; (void) w28;
    const double w35 = P.w35; (void) w35;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    (void) t; (void) q;
    const double w0 = std::pow(q, 2);
    const double w1 = ik*w0;
    const double w2 = 6*q;
    const double w3 = ik*w2;
    const double w4 = w0 - w3;
    const double w5 = ik - w2;
    const double w6 = w1 + w4 + w5 + 1;
    const double w8 = std::pow(t, 4);
    const double w9 = w8/std::pow(s, 4);
    const double w11 = ik*q;
    const double w12 = q + w10 + w11 - 1;
    const double w13 = 24*w9*z;
    const double w16 = 2*w8;
    const double w17 = w15*w16;
    const double w18 = 6*w9;
    const double w19 = 6*w0;
    const double w20 = -8*ik*q;
    const double w21 = ik*w17;
    const double w25 = 2*q;
    const double w26 = 2*w11;
    const double w27 = -w0;
    const double w29 = w28*w8;
    const double w30 = (-ik + q)/s;
    const double w31 = std::pow(z, 2);
    const double w32 = 5*q;
    const double w33 = std::pow(q, 3);
    const double w34 = 4*w0;
    const double w36 = w35*w8;
    const double w37 = 12*w36;
    const double w38 = std::pow(q, 4)*w8;
    o_mu_mu_mu_mu[i] = w6*w7*w9;
    o_mu_mu_mu_sigma[i] = -ik*w12*w13;
    o_mu_mu_mu_nu[i] = w14*w17*z*(-8*q + w4 + w7 + 3);
    o_mu_mu_sigma_sigma[i] = -w18*w6;
    o_mu_mu_sigma_nu[i] = w21*(-3*w1 - w19 - w20 - w5);
    o_mu_mu_nu_nu[i] = w16*w23*(3*q - 3*w0 + w10 + 5*w11);
    o_mu_sigma_sigma_sigma[i] = w12*w13;
    o_mu_sigma_sigma_nu[i] = w21*z*(10*q - 2*w0 - w1 - w20 - w24);
    o_mu_sigma_nu_nu[i] = 4*w23*w8*z*(w10 + w25 + w26 + w27);
    o_mu_nu_nu_nu[i] = 6*w29*w30*z;
    o_sigma_sigma_sigma_sigma[i] = w18*(-w2 - w27 - w31*w32 - w31*w33 - w31*w34 - 10*w31 + 1);
    o_sigma_sigma_sigma_nu[i] = w15*w25*w8*(15*q + w3 + w33 + w34 - w7);
    o_sigma_sigma_nu_nu[i] = w22*w26*w8*(3*w11 - w24 + w27 + w32);
    o_sigma_nu_nu_nu[i] = w14*w2*w30*w8;
    o_nu_nu_nu_nu[i] = V5 - q*w37 - 18*w0*w36 + w19*w29 + w28*w38 + 4*w29*w33 - w33*w37 - 3*w35*w38;
  });
  return List::create(Named("mu_mu_mu_mu") = o_mu_mu_mu_mu, Named("mu_mu_mu_sigma") = o_mu_mu_mu_sigma, Named("mu_mu_mu_nu") = o_mu_mu_mu_nu, Named("mu_mu_sigma_sigma") = o_mu_mu_sigma_sigma, Named("mu_mu_sigma_nu") = o_mu_mu_sigma_nu, Named("mu_mu_nu_nu") = o_mu_mu_nu_nu, Named("mu_sigma_sigma_sigma") = o_mu_sigma_sigma_sigma, Named("mu_sigma_sigma_nu") = o_mu_sigma_sigma_nu, Named("mu_sigma_nu_nu") = o_mu_sigma_nu_nu, Named("mu_nu_nu_nu") = o_mu_nu_nu_nu, Named("sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma, Named("sigma_sigma_sigma_nu") = o_sigma_sigma_sigma_nu, Named("sigma_sigma_nu_nu") = o_sigma_sigma_nu_nu, Named("sigma_nu_nu_nu") = o_sigma_nu_nu_nu, Named("nu_nu_nu_nu") = o_nu_nu_nu_nu);
}

// order 4, expected
// [[Rcpp::export]]
List student_t1_deriv4_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu_mu_mu(n);
  NumericVector o_mu_mu_mu_sigma(n);
  NumericVector o_mu_mu_mu_nu(n);
  NumericVector o_mu_mu_sigma_sigma(n);
  NumericVector o_mu_mu_sigma_nu(n);
  NumericVector o_mu_mu_nu_nu(n);
  NumericVector o_mu_sigma_sigma_sigma(n);
  NumericVector o_mu_sigma_sigma_nu(n);
  NumericVector o_mu_sigma_nu_nu(n);
  NumericVector o_mu_nu_nu_nu(n);
  NumericVector o_sigma_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_sigma_nu(n);
  NumericVector o_sigma_sigma_nu_nu(n);
  NumericVector o_sigma_nu_nu_nu(n);
  NumericVector o_nu_nu_nu_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, V6, w0, w1, w2, w3, w4, w5, w6, w7, w8, w9; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double V6 = t1_V6(v);
    const double w0 = std::pow(ik, 2);
    const double w1 = 2*w0;
    const double w2 = 1/(std::pow(s, 4)*(12*ik + 35*w0 + 1));
    const double w3 = 6*w2*(3*ik + w1 + 1);
    const double w4 = std::pow(ik, 3);
    const double w5 = 1.0/(15*ik + 71*w0 + 105*w4 + 1);
    const double w6 = w5/std::pow(s, 3);
    const double w7 = 4*w0 - 1;
    const double w8 = w4*w5/std::pow(s, 2);
    const double w9 = std::pow(ik, 4);
    return Prm{m, s, ik, V6, w0, w1, w2, w3, w4, w5, w6, w7, w8, w9};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double V6 = P.V6; (void) V6;
    const double w0 = P.w0; (void) w0;
    const double w1 = P.w1; (void) w1;
    const double w2 = P.w2; (void) w2;
    const double w3 = P.w3; (void) w3;
    const double w4 = P.w4; (void) w4;
    const double w5 = P.w5; (void) w5;
    const double w6 = P.w6; (void) w6;
    const double w7 = P.w7; (void) w7;
    const double w8 = P.w8; (void) w8;
    const double w9 = P.w9; (void) w9;
    o_mu_mu_mu_mu[i] = ik*w3;
    o_mu_mu_mu_sigma[i] = 0;
    o_mu_mu_mu_nu[i] = 0;
    o_mu_mu_sigma_sigma[i] = -w3;
    o_mu_mu_sigma_nu[i] = w1*w6*(9*ik - w1 + 5);
    o_mu_mu_nu_nu[i] = -4*w7*w8;
    o_mu_sigma_sigma_sigma[i] = 0;
    o_mu_sigma_sigma_nu[i] = 0;
    o_mu_sigma_nu_nu[i] = 0;
    o_mu_nu_nu_nu[i] = 0;
    o_sigma_sigma_sigma_sigma[i] = -18*w2*(11*ik + 3);
    o_sigma_sigma_sigma_nu[i] = 6*w0*w6*(31*ik + 13);
    o_sigma_sigma_nu_nu[i] = -12*w8*(ik - 2);
    o_sigma_nu_nu_nu[i] = -12*w7*w9/(s*(16*ik + 86*w0 + 176*w4 + 105*w9 + 1));
    o_nu_nu_nu_nu[i] = V6;
  });
  return List::create(Named("mu_mu_mu_mu") = o_mu_mu_mu_mu, Named("mu_mu_mu_sigma") = o_mu_mu_mu_sigma, Named("mu_mu_mu_nu") = o_mu_mu_mu_nu, Named("mu_mu_sigma_sigma") = o_mu_mu_sigma_sigma, Named("mu_mu_sigma_nu") = o_mu_mu_sigma_nu, Named("mu_mu_nu_nu") = o_mu_mu_nu_nu, Named("mu_sigma_sigma_sigma") = o_mu_sigma_sigma_sigma, Named("mu_sigma_sigma_nu") = o_mu_sigma_sigma_nu, Named("mu_sigma_nu_nu") = o_mu_sigma_nu_nu, Named("mu_nu_nu_nu") = o_mu_nu_nu_nu, Named("sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma, Named("sigma_sigma_sigma_nu") = o_sigma_sigma_sigma_nu, Named("sigma_sigma_nu_nu") = o_sigma_sigma_nu_nu, Named("sigma_nu_nu_nu") = o_sigma_nu_nu_nu, Named("nu_nu_nu_nu") = o_nu_nu_nu_nu);
}

// order 5 of log f in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t1_deriv5_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu_mu_mu_mu(n);
  NumericVector o_mu_mu_mu_mu_sigma(n);
  NumericVector o_mu_mu_mu_mu_nu(n);
  NumericVector o_mu_mu_mu_sigma_sigma(n);
  NumericVector o_mu_mu_mu_sigma_nu(n);
  NumericVector o_mu_mu_mu_nu_nu(n);
  NumericVector o_mu_mu_sigma_sigma_sigma(n);
  NumericVector o_mu_mu_sigma_sigma_nu(n);
  NumericVector o_mu_mu_sigma_nu_nu(n);
  NumericVector o_mu_mu_nu_nu_nu(n);
  NumericVector o_mu_sigma_sigma_sigma_sigma(n);
  NumericVector o_mu_sigma_sigma_sigma_nu(n);
  NumericVector o_mu_sigma_sigma_nu_nu(n);
  NumericVector o_mu_sigma_nu_nu_nu(n);
  NumericVector o_mu_nu_nu_nu_nu(n);
  NumericVector o_sigma_sigma_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_sigma_sigma_nu(n);
  NumericVector o_sigma_sigma_sigma_nu_nu(n);
  NumericVector o_sigma_sigma_nu_nu_nu(n);
  NumericVector o_sigma_nu_nu_nu_nu(n);
  NumericVector o_nu_nu_nu_nu_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, V7, w0, w3, w15, w22, w23, w32, w33, w34, w37, w42, w45; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double V7 = t1_V7(v);
    const double w0 = std::pow(ik, 2);
    const double w3 = 5*ik;
    const double w15 = 2*ik;
    const double w22 = std::pow(s, -4);
    const double w23 = 6*w22;
    const double w32 = std::pow(s, -3);
    const double w33 = 3*ik;
    const double w34 = std::pow(ik, 3);
    const double w37 = -ik;
    const double w42 = std::pow(s, -2);
    const double w45 = std::pow(ik, 4);
    return Prm{m, s, ik, V7, w0, w3, w15, w22, w23, w32, w33, w34, w37, w42, w45};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double V7 = P.V7; (void) V7;
    const double w0 = P.w0; (void) w0;
    const double w3 = P.w3; (void) w3;
    const double w15 = P.w15; (void) w15;
    const double w22 = P.w22; (void) w22;
    const double w23 = P.w23; (void) w23;
    const double w32 = P.w32; (void) w32;
    const double w33 = P.w33; (void) w33;
    const double w34 = P.w34; (void) w34;
    const double w37 = P.w37; (void) w37;
    const double w42 = P.w42; (void) w42;
    const double w45 = P.w45; (void) w45;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    (void) t; (void) q;
    const double w1 = std::pow(q, 2);
    const double w2 = ik*w1;
    const double w4 = 10*q;
    const double w5 = -w4;
    const double w6 = -ik*w4 + w5;
    const double w7 = w1 + w2 + w3 + w6 + 5;
    const double w8 = std::pow(t, 5);
    const double w9 = 24*w8;
    const double w10 = w9/std::pow(s, 5);
    const double w11 = w10*w7*z;
    const double w12 = 5*w1;
    const double w13 = ik + w1*w3 + w12 + w6 + 1;
    const double w14 = ik*w10;
    const double w16 = std::pow(q, 3);
    const double w17 = -w16;
    const double w18 = ik*q;
    const double w19 = 20*w18;
    const double w20 = 15*w1 + 10*w2;
    const double w21 = w0*w8;
    const double w24 = 2*w1;
    const double w25 = 5*q;
    const double w26 = q*w3;
    const double w27 = -w25 - w26;
    const double w28 = w24 + w27 + 1;
    const double w29 = w15 + w2;
    const double w30 = w9*z;
    const double w31 = w22*w30;
    const double w35 = w34*w8;
    const double w36 = 12*z;
    const double w38 = 15*ik*q;
    const double w39 = 6*q;
    const double w40 = 4*w21*w32;
    const double w41 = -w12;
    const double w43 = w35*w42;
    const double w44 = 3*w1 + w27;
    const double w46 = (ik - q)/s;
    const double w47 = std::pow(z, 2);
    const double w48 = std::pow(q, 4);
    const double w49 = 5*w16;
    const double w50 = std::pow(ik, 5)*w8;
    const double w51 = 60*w50;
    const double w52 = w45*w8;
    const double w53 = 30*w52;
    const double w54 = std::pow(q, 5);
    const double w55 = 120*w50;
    o_mu_mu_mu_mu_mu[i] = w0*w11;
    o_mu_mu_mu_mu_sigma[i] = -w13*w14;
    o_mu_mu_mu_mu_nu[i] = w21*w23*(15*q - w15 - w17 + w19 - w20 - 1);
    o_mu_mu_mu_sigma_sigma[i] = -w14*w7*z;
    o_mu_mu_mu_sigma_nu[i] = w0*w31*(-w28 - w29);
    o_mu_mu_mu_nu_nu[i] = w32*w35*w36*(-w28 - w33);
    o_mu_mu_sigma_sigma_sigma[i] = w10*w13;
    o_mu_mu_sigma_sigma_nu[i] = ik*w23*w8*(15*ik*w1 - ik*w16 + 20*w1 - 2*w16 - w37 - w38 - w4);
    o_mu_mu_sigma_nu_nu[i] = w40*(ik - 3*w16 - 13*w18 + w20 - w39);
    o_mu_mu_nu_nu_nu[i] = 6*w43*(-3*q - 7*w18 - w37 - w41);
    o_mu_sigma_sigma_sigma_sigma[i] = w11;
    o_mu_sigma_sigma_sigma_nu[i] = ik*w31*(ik + 2*w2 + w44);
    o_mu_sigma_sigma_nu_nu[i] = w40*z*(13*w1 + w17 + 6*w2 + w33 - w38 + w5);
    o_mu_sigma_nu_nu_nu[i] = w36*w43*(-q*w33 - 2*q + w24 - w37);
    o_mu_nu_nu_nu_nu[i] = w30*w45*w46;
    o_sigma_sigma_sigma_sigma_sigma[i] = w10*(10*w1*w47 + w25*w47 + w4 + w41 + w47*w48 + w47*w49 + 15*w47 - 1);
    o_sigma_sigma_sigma_sigma_nu[i] = w22*w39*w8*(2*ik*w1 + 10*ik - 35*q - 7*w1 - w19 - w48 - w49);
    o_sigma_sigma_sigma_nu_nu[i] = 12*w18*w32*w8*(w29 + w44);
    o_sigma_sigma_nu_nu_nu[i] = w21*w39*w42*(3*w1 - w25 - w26 + w33);
    o_sigma_nu_nu_nu_nu[i] = q*w34*w46*w9;
    o_nu_nu_nu_nu_nu[i] = V7 + q*w51 - w1*w53 + w1*w55 - w16*w53 + w16*w55 + w48*w51 - 15*w48*w52 + 12*w50*w54 - 3*w52*w54;
  });
  return List::create(Named("mu_mu_mu_mu_mu") = o_mu_mu_mu_mu_mu, Named("mu_mu_mu_mu_sigma") = o_mu_mu_mu_mu_sigma, Named("mu_mu_mu_mu_nu") = o_mu_mu_mu_mu_nu, Named("mu_mu_mu_sigma_sigma") = o_mu_mu_mu_sigma_sigma, Named("mu_mu_mu_sigma_nu") = o_mu_mu_mu_sigma_nu, Named("mu_mu_mu_nu_nu") = o_mu_mu_mu_nu_nu, Named("mu_mu_sigma_sigma_sigma") = o_mu_mu_sigma_sigma_sigma, Named("mu_mu_sigma_sigma_nu") = o_mu_mu_sigma_sigma_nu, Named("mu_mu_sigma_nu_nu") = o_mu_mu_sigma_nu_nu, Named("mu_mu_nu_nu_nu") = o_mu_mu_nu_nu_nu, Named("mu_sigma_sigma_sigma_sigma") = o_mu_sigma_sigma_sigma_sigma, Named("mu_sigma_sigma_sigma_nu") = o_mu_sigma_sigma_sigma_nu, Named("mu_sigma_sigma_nu_nu") = o_mu_sigma_sigma_nu_nu, Named("mu_sigma_nu_nu_nu") = o_mu_sigma_nu_nu_nu, Named("mu_nu_nu_nu_nu") = o_mu_nu_nu_nu_nu, Named("sigma_sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma_sigma, Named("sigma_sigma_sigma_sigma_nu") = o_sigma_sigma_sigma_sigma_nu, Named("sigma_sigma_sigma_nu_nu") = o_sigma_sigma_sigma_nu_nu, Named("sigma_sigma_nu_nu_nu") = o_sigma_sigma_nu_nu_nu, Named("sigma_nu_nu_nu_nu") = o_sigma_nu_nu_nu_nu, Named("nu_nu_nu_nu_nu") = o_nu_nu_nu_nu_nu);
}

// first derivatives of the expected Hessian
// [[Rcpp::export]]
List student_t1_dexpected1_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu_mu(n);
  NumericVector o_mu_mu_sigma(n);
  NumericVector o_mu_mu_nu(n);
  NumericVector o_sigma_sigma_mu(n);
  NumericVector o_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_nu(n);
  NumericVector o_nu_nu_mu(n);
  NumericVector o_nu_nu_sigma(n);
  NumericVector o_nu_nu_nu(n);
  NumericVector o_mu_sigma_mu(n);
  NumericVector o_mu_sigma_sigma(n);
  NumericVector o_mu_sigma_nu(n);
  NumericVector o_mu_nu_mu(n);
  NumericVector o_mu_nu_sigma(n);
  NumericVector o_mu_nu_nu(n);
  NumericVector o_sigma_nu_mu(n);
  NumericVector o_sigma_nu_sigma(n);
  NumericVector o_sigma_nu_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, V8, w0, w1, w2, w3, w4; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double V8 = t1_V8(v);
    const double w0 = 1/(std::pow(s, 3)*(3*ik + 1));
    const double w1 = std::pow(ik, 2);
    const double w2 = std::pow(s, -2);
    const double w3 = w1*w2/(6*ik + 9*w1 + 1);
    const double w4 = std::pow(ik, 3);
    return Prm{m, s, ik, V8, w0, w1, w2, w3, w4};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double V8 = P.V8; (void) V8;
    const double w0 = P.w0; (void) w0;
    const double w1 = P.w1; (void) w1;
    const double w2 = P.w2; (void) w2;
    const double w3 = P.w3; (void) w3;
    const double w4 = P.w4; (void) w4;
    o_mu_mu_mu[i] = 0;
    o_mu_mu_sigma[i] = 2*w0*(ik + 1);
    o_mu_mu_nu[i] = -2*w3;
    o_sigma_sigma_mu[i] = 0;
    o_sigma_sigma_sigma[i] = 4*w0;
    o_sigma_sigma_nu[i] = -6*w3;
    o_nu_nu_mu[i] = 0;
    o_nu_nu_sigma[i] = 0;
    o_nu_nu_nu[i] = V8;
    o_mu_sigma_mu[i] = 0;
    o_mu_sigma_sigma[i] = 0;
    o_mu_sigma_nu[i] = 0;
    o_mu_nu_mu[i] = 0;
    o_mu_nu_sigma[i] = 0;
    o_mu_nu_nu[i] = 0;
    o_sigma_nu_mu[i] = 0;
    o_sigma_nu_sigma[i] = -2*w1*w2/(4*ik + 3*w1 + 1);
    o_sigma_nu_nu[i] = -4*w4*(2*ik + 1)/(s*(9*std::pow(ik, 4) + 8*ik + 22*w1 + 24*w4 + 1));
  });
  return List::create(Named("mu_mu_mu") = o_mu_mu_mu, Named("mu_mu_sigma") = o_mu_mu_sigma, Named("mu_mu_nu") = o_mu_mu_nu, Named("sigma_sigma_mu") = o_sigma_sigma_mu, Named("sigma_sigma_sigma") = o_sigma_sigma_sigma, Named("sigma_sigma_nu") = o_sigma_sigma_nu, Named("nu_nu_mu") = o_nu_nu_mu, Named("nu_nu_sigma") = o_nu_nu_sigma, Named("nu_nu_nu") = o_nu_nu_nu, Named("mu_sigma_mu") = o_mu_sigma_mu, Named("mu_sigma_sigma") = o_mu_sigma_sigma, Named("mu_sigma_nu") = o_mu_sigma_nu, Named("mu_nu_mu") = o_mu_nu_mu, Named("mu_nu_sigma") = o_mu_nu_sigma, Named("mu_nu_nu") = o_mu_nu_nu, Named("sigma_nu_mu") = o_sigma_nu_mu, Named("sigma_nu_sigma") = o_sigma_nu_sigma, Named("sigma_nu_nu") = o_sigma_nu_nu);
}

// second derivatives of the expected Hessian
// [[Rcpp::export]]
List student_t1_dexpected2_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu_mu_mu(n);
  NumericVector o_mu_mu_sigma_sigma(n);
  NumericVector o_mu_mu_nu_nu(n);
  NumericVector o_mu_mu_mu_sigma(n);
  NumericVector o_mu_mu_mu_nu(n);
  NumericVector o_mu_mu_sigma_nu(n);
  NumericVector o_sigma_sigma_mu_mu(n);
  NumericVector o_sigma_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_nu_nu(n);
  NumericVector o_sigma_sigma_mu_sigma(n);
  NumericVector o_sigma_sigma_mu_nu(n);
  NumericVector o_sigma_sigma_sigma_nu(n);
  NumericVector o_nu_nu_mu_mu(n);
  NumericVector o_nu_nu_sigma_sigma(n);
  NumericVector o_nu_nu_nu_nu(n);
  NumericVector o_nu_nu_mu_sigma(n);
  NumericVector o_nu_nu_mu_nu(n);
  NumericVector o_nu_nu_sigma_nu(n);
  NumericVector o_mu_sigma_mu_mu(n);
  NumericVector o_mu_sigma_sigma_sigma(n);
  NumericVector o_mu_sigma_nu_nu(n);
  NumericVector o_mu_sigma_mu_sigma(n);
  NumericVector o_mu_sigma_mu_nu(n);
  NumericVector o_mu_sigma_sigma_nu(n);
  NumericVector o_mu_nu_mu_mu(n);
  NumericVector o_mu_nu_sigma_sigma(n);
  NumericVector o_mu_nu_nu_nu(n);
  NumericVector o_mu_nu_mu_sigma(n);
  NumericVector o_mu_nu_mu_nu(n);
  NumericVector o_mu_nu_sigma_nu(n);
  NumericVector o_sigma_nu_mu_mu(n);
  NumericVector o_sigma_nu_sigma_sigma(n);
  NumericVector o_sigma_nu_nu_nu(n);
  NumericVector o_sigma_nu_mu_sigma(n);
  NumericVector o_sigma_nu_mu_nu(n);
  NumericVector o_sigma_nu_sigma_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, V9, w0, w1, w2, w3, w4, w5, w6, w7, w8; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double V9 = t1_V9(v);
    const double w0 = 1/(std::pow(s, 4)*(3*ik + 1));
    const double w1 = std::pow(ik, 3);
    const double w2 = std::pow(s, -2);
    const double w3 = std::pow(ik, 2);
    const double w4 = w1*w2/(9*ik + 27*w1 + 27*w3 + 1);
    const double w5 = std::pow(s, -3);
    const double w6 = w3*w5/(6*ik + 9*w3 + 1);
    const double w7 = std::pow(ik, 4);
    const double w8 = 12*ik;
    return Prm{m, s, ik, V9, w0, w1, w2, w3, w4, w5, w6, w7, w8};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double V9 = P.V9; (void) V9;
    const double w0 = P.w0; (void) w0;
    const double w1 = P.w1; (void) w1;
    const double w2 = P.w2; (void) w2;
    const double w3 = P.w3; (void) w3;
    const double w4 = P.w4; (void) w4;
    const double w5 = P.w5; (void) w5;
    const double w6 = P.w6; (void) w6;
    const double w7 = P.w7; (void) w7;
    const double w8 = P.w8; (void) w8;
    o_mu_mu_mu_mu[i] = 0;
    o_mu_mu_sigma_sigma[i] = -6*w0*(ik + 1);
    o_mu_mu_nu_nu[i] = 4*w4;
    o_mu_mu_mu_sigma[i] = 0;
    o_mu_mu_mu_nu[i] = 0;
    o_mu_mu_sigma_nu[i] = 4*w6;
    o_sigma_sigma_mu_mu[i] = 0;
    o_sigma_sigma_sigma_sigma[i] = -12*w0;
    o_sigma_sigma_nu_nu[i] = 12*w4;
    o_sigma_sigma_mu_sigma[i] = 0;
    o_sigma_sigma_mu_nu[i] = 0;
    o_sigma_sigma_sigma_nu[i] = 12*w6;
    o_nu_nu_mu_mu[i] = 0;
    o_nu_nu_sigma_sigma[i] = 0;
    o_nu_nu_nu_nu[i] = V9;
    o_nu_nu_mu_sigma[i] = 0;
    o_nu_nu_mu_nu[i] = 0;
    o_nu_nu_sigma_nu[i] = 0;
    o_mu_sigma_mu_mu[i] = 0;
    o_mu_sigma_sigma_sigma[i] = 0;
    o_mu_sigma_nu_nu[i] = 0;
    o_mu_sigma_mu_sigma[i] = 0;
    o_mu_sigma_mu_nu[i] = 0;
    o_mu_sigma_sigma_nu[i] = 0;
    o_mu_nu_mu_mu[i] = 0;
    o_mu_nu_sigma_sigma[i] = 0;
    o_mu_nu_nu_nu[i] = 0;
    o_mu_nu_mu_sigma[i] = 0;
    o_mu_nu_mu_nu[i] = 0;
    o_mu_nu_sigma_nu[i] = 0;
    o_sigma_nu_mu_mu[i] = 0;
    o_sigma_nu_sigma_sigma[i] = 4*w3*w5/(4*ik + 3*w3 + 1);
    o_sigma_nu_nu_nu[i] = 4*w7*(13*w3 + w8 + 3)/(s*(27*std::pow(ik, 6) + 108*std::pow(ik, 5) + 136*w1 + 57*w3 + 171*w7 + w8 + 1));
    o_sigma_nu_mu_sigma[i] = 0;
    o_sigma_nu_mu_nu[i] = 0;
    o_sigma_nu_sigma_nu[i] = 4*w1*w2*(2*ik + 1)/(8*ik + 24*w1 + 22*w3 + 9*w7 + 1);
  });
  return List::create(Named("mu_mu_mu_mu") = o_mu_mu_mu_mu, Named("mu_mu_sigma_sigma") = o_mu_mu_sigma_sigma, Named("mu_mu_nu_nu") = o_mu_mu_nu_nu, Named("mu_mu_mu_sigma") = o_mu_mu_mu_sigma, Named("mu_mu_mu_nu") = o_mu_mu_mu_nu, Named("mu_mu_sigma_nu") = o_mu_mu_sigma_nu, Named("sigma_sigma_mu_mu") = o_sigma_sigma_mu_mu, Named("sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma, Named("sigma_sigma_nu_nu") = o_sigma_sigma_nu_nu, Named("sigma_sigma_mu_sigma") = o_sigma_sigma_mu_sigma, Named("sigma_sigma_mu_nu") = o_sigma_sigma_mu_nu, Named("sigma_sigma_sigma_nu") = o_sigma_sigma_sigma_nu, Named("nu_nu_mu_mu") = o_nu_nu_mu_mu, Named("nu_nu_sigma_sigma") = o_nu_nu_sigma_sigma, Named("nu_nu_nu_nu") = o_nu_nu_nu_nu, Named("nu_nu_mu_sigma") = o_nu_nu_mu_sigma, Named("nu_nu_mu_nu") = o_nu_nu_mu_nu, Named("nu_nu_sigma_nu") = o_nu_nu_sigma_nu, Named("mu_sigma_mu_mu") = o_mu_sigma_mu_mu, Named("mu_sigma_sigma_sigma") = o_mu_sigma_sigma_sigma, Named("mu_sigma_nu_nu") = o_mu_sigma_nu_nu, Named("mu_sigma_mu_sigma") = o_mu_sigma_mu_sigma, Named("mu_sigma_mu_nu") = o_mu_sigma_mu_nu, Named("mu_sigma_sigma_nu") = o_mu_sigma_sigma_nu, Named("mu_nu_mu_mu") = o_mu_nu_mu_mu, Named("mu_nu_sigma_sigma") = o_mu_nu_sigma_sigma, Named("mu_nu_nu_nu") = o_mu_nu_nu_nu, Named("mu_nu_mu_sigma") = o_mu_nu_mu_sigma, Named("mu_nu_mu_nu") = o_mu_nu_mu_nu, Named("mu_nu_sigma_nu") = o_mu_nu_sigma_nu, Named("sigma_nu_mu_mu") = o_sigma_nu_mu_mu, Named("sigma_nu_sigma_sigma") = o_sigma_nu_sigma_sigma, Named("sigma_nu_nu_nu") = o_sigma_nu_nu_nu, Named("sigma_nu_mu_sigma") = o_sigma_nu_mu_sigma, Named("sigma_nu_mu_nu") = o_sigma_nu_mu_nu, Named("sigma_nu_sigma_nu") = o_sigma_nu_sigma_nu);
}

// 1 in y and 1 in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t1_cross_y_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu(n);
  NumericVector o_sigma(n);
  NumericVector o_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    return Prm{m, s, ik};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    (void) t; (void) q;
    const double w0 = std::pow(t, 2);
    const double w1 = w0/std::pow(s, 2);
    o_mu[i] = w1*(-ik*q + ik - q + 1);
    o_sigma[i] = 2*w1*z*(ik + 1);
    o_nu[i] = ik*w0*z*(ik - q)/s;
  });
  return List::create(Named("mu") = o_mu, Named("sigma") = o_sigma, Named("nu") = o_nu);
}

// 2 in y and 1 in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t1_cross2_y_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu(n);
  NumericVector o_sigma(n);
  NumericVector o_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, w0; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double w0 = 3*ik;
    return Prm{m, s, ik, w0};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double w0 = P.w0; (void) w0;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    (void) t; (void) q;
    const double w1 = std::pow(t, 3);
    const double w2 = 2*w1/std::pow(s, 3);
    const double w3 = -ik + q*w0 + 3*q;
    o_mu[i] = ik*w2*z*(ik*q + q - w0 - 3);
    o_sigma[i] = w2*(1 - w3);
    o_nu[i] = ik*w1*(std::pow(q, 2) - w3)/std::pow(s, 2);
  });
  return List::create(Named("mu") = o_mu, Named("sigma") = o_sigma, Named("nu") = o_nu);
}

// 1 in y and 2 in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t1_grad_y_hess_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu(n);
  NumericVector o_sigma_sigma(n);
  NumericVector o_nu_nu(n);
  NumericVector o_mu_sigma(n);
  NumericVector o_mu_nu(n);
  NumericVector o_sigma_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, w7, w9; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double w7 = -ik;
    const double w9 = ik/std::pow(s, 2);
    return Prm{m, s, ik, w7, w9};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double w7 = P.w7; (void) w7;
    const double w9 = P.w9; (void) w9;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    (void) t; (void) q;
    const double w0 = ik*q;
    const double w1 = -3*ik + q + w0 - 3;
    const double w2 = std::pow(t, 3);
    const double w3 = 2*w2;
    const double w4 = w3/std::pow(s, 3);
    const double w5 = w4*z;
    const double w6 = w3*z;
    const double w8 = 3*ik*q + 3*q + w7;
    o_mu_mu[i] = -ik*w1*w5;
    o_sigma_sigma[i] = w1*w5;
    o_nu_nu[i] = std::pow(ik, 2)*w6*(-ik + q)/s;
    o_mu_sigma[i] = w4*(w8 - 1);
    o_mu_nu[i] = w2*w9*(-std::pow(q, 2) + w8);
    o_sigma_nu[i] = w6*w9*(2*q + w0 + w7);
  });
  return List::create(Named("mu_mu") = o_mu_mu, Named("sigma_sigma") = o_sigma_sigma, Named("nu_nu") = o_nu_nu, Named("mu_sigma") = o_mu_sigma, Named("mu_nu") = o_mu_nu, Named("sigma_nu") = o_sigma_nu);
}

// 2 in y and 2 in (mu, sigma, nu)
// [[Rcpp::export]]
List student_t1_hess_y_hess_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_mu_mu(n);
  NumericVector o_sigma_sigma(n);
  NumericVector o_nu_nu(n);
  NumericVector o_mu_sigma(n);
  NumericVector o_mu_nu(n);
  NumericVector o_sigma_nu(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik, w6, w9, w13; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    const double w6 = 6*ik;
    const double w9 = -ik;
    const double w13 = std::pow(s, -3);
    return Prm{m, s, ik, w6, w9, w13};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double w6 = P.w6; (void) w6;
    const double w9 = P.w9; (void) w9;
    const double w13 = P.w13; (void) w13;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    (void) t; (void) q;
    const double w0 = std::pow(q, 2);
    const double w1 = ik*w0;
    const double w2 = 6*q;
    const double w3 = -ik*w2 + w0;
    const double w4 = ik - w2;
    const double w5 = w1 + w3 + w4 + 1;
    const double w7 = std::pow(t, 4);
    const double w8 = w7/std::pow(s, 4);
    const double w10 = ik*q;
    const double w11 = 2*w7;
    const double w12 = std::pow(ik, 2)*w11;
    o_mu_mu[i] = w5*w6*w8;
    o_sigma_sigma[i] = -6*w5*w8;
    o_nu_nu[i] = w12*(3*q - 3*w0 + 5*w10 + w9)/std::pow(s, 2);
    o_mu_sigma[i] = 24*ik*w8*z*(-q - w10 - w9 + 1);
    o_mu_nu[i] = w12*w13*z*(-8*q + w3 + w6 + 3);
    o_sigma_nu[i] = ik*w11*w13*(8*ik*q - 6*w0 - 3*w1 - w4);
  });
  return List::create(Named("mu_mu") = o_mu_mu, Named("sigma_sigma") = o_sigma_sigma, Named("nu_nu") = o_nu_nu, Named("mu_sigma") = o_mu_sigma, Named("mu_nu") = o_mu_nu, Named("sigma_nu") = o_sigma_nu);
}

// order 1 of log f in y
// [[Rcpp::export]]
List student_t1_dy1_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_y(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    return Prm{m, s, ik};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    (void) t; (void) q;
    o_y[i] = -t*z*(ik + 1)/s;
  });
  return List::create(Named("y") = o_y);
}

// order 2 of log f in y
// [[Rcpp::export]]
List student_t1_dy2_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_y(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    return Prm{m, s, ik};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    (void) t; (void) q;
    o_y[i] = std::pow(t, 2)*(ik*q - ik + q - 1)/std::pow(s, 2);
  });
  return List::create(Named("y") = o_y);
}

// order 3 of log f in y
// [[Rcpp::export]]
List student_t1_dy3_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_y(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    return Prm{m, s, ik};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    (void) t; (void) q;
    o_y[i] = 2*ik*std::pow(t, 3)*z*(-ik*q + 3*ik - q + 3)/std::pow(s, 3);
  });
  return List::create(Named("y") = o_y);
}

// order 4 of log f in y
// [[Rcpp::export]]
List student_t1_dy4_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu, int threads = 1) {
  const int n = y.size();
  const int n_mu = mu.size(), n_sigma = sigma.size(), n_nu = nu.size();
  NumericVector o_y(n);
  // the quantities of the parameters alone, once when they are scalars
  struct Prm { double m, s, ik; };
  auto make = [&](std::size_t i) {
    const double m = mu[i % n_mu], s = sigma[i % n_sigma];
    const double v = nu[i % n_nu];
    const double ik = 1.0 / (v - 0.0);
    return Prm{m, s, ik};
  };
  const bool scalar = n_mu == 1 && n_sigma == 1 && n_nu == 1;
  Prm P0{};
  if (scalar) P0 = make(0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const Prm P = scalar ? P0 : make(i);
    const double m = P.m; (void) m;
    const double s = P.s; (void) s;
    const double ik = P.ik; (void) ik;
    const double z = (y[i] - m) / s;
    const double q = z * z * ik;
    const double t = 1.0 / (1.0 + q);
    (void) t; (void) q;
    const double w0 = 6*q;
    const double w1 = std::pow(q, 2);
    o_y[i] = 6*ik*std::pow(t, 4)*(-ik*w0 + ik*w1 + ik - w0 + w1 + 1)/std::pow(s, 4);
  });
  return List::create(Named("y") = o_y);
}
