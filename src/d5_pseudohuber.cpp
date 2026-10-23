#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "pt_pseudohuber.h"
#include "pt_ph_nu.h"
using namespace Rcpp;

// The fifth derivatives of the pseudo-Huber log-density in (mu, sigma, nu):
//   l = -D - log 2 - log sigma - B(nu),  D = sqrt(nu + z^2), z = (y - mu)/sigma,
// with B(nu) = log(sqrt(nu) K1(sqrt(nu))) the nu-only part of pt_ph_nu.h.
// The components are written in p = z/D (bounded by one) and D, as the
// fourth order is written in rho and kappa^2 (pseudohuber_hd.cpp), so no
// power of z is formed; the pure-nu component adds B^(5)(nu) from
// ph_nu_derivs(), to order five. The loop body is derived in
// misure_gas_congiunto/c/fase4/D5/pseudohuber/derive.py.

// [[Rcpp::export]]
List pseudohuber_deriv5_cpp(NumericVector y, NumericVector mu,
                            NumericVector sigma, NumericVector nu) {
  const int n = std::max(y.size(), std::max(mu.size(), std::max(sigma.size(), nu.size())));
  const int n_y = y.size(), n_mu = mu.size(), n_s = sigma.size(), n_v = nu.size();
  NumericVector o_mu_mu_mu_mu_mu(n), o_mu_mu_mu_mu_sigma(n), o_mu_mu_mu_mu_nu(n),
    o_mu_mu_mu_sigma_sigma(n), o_mu_mu_mu_sigma_nu(n), o_mu_mu_mu_nu_nu(n),
    o_mu_mu_sigma_sigma_sigma(n), o_mu_mu_sigma_sigma_nu(n), o_mu_mu_sigma_nu_nu(n),
    o_mu_mu_nu_nu_nu(n), o_mu_sigma_sigma_sigma_sigma(n), o_mu_sigma_sigma_sigma_nu(n),
    o_mu_sigma_sigma_nu_nu(n), o_mu_sigma_nu_nu_nu(n), o_mu_nu_nu_nu_nu(n),
    o_sigma_sigma_sigma_sigma_sigma(n), o_sigma_sigma_sigma_sigma_nu(n),
    o_sigma_sigma_sigma_nu_nu(n), o_sigma_sigma_nu_nu_nu(n), o_sigma_nu_nu_nu_nu(n),
    o_nu_nu_nu_nu_nu(n);
  for (int i = 0; i < n; i++) {
    const double yy = y[i % n_y];
    const double mm = mu[i % n_mu];
    const double ss = sigma[i % n_s];
    const double nn = nu[i % n_v];
    const double zz = (yy - mm)/ss;
    const double DD = std::hypot(std::sqrt(nn), zz);
    const double pp = zz/DD;
    double dB[6];
    double dlR[6];
    double R0;
    d7::ph_nu_derivs(nn, 5, dB, dlR, R0);
    const double B5 = dB[5];
        const double t0 = std::pow(pp, 2);
        const double t1 = 7*t0;
        const double t2 = t1 - 3;
        const double t3 = pp/std::pow(DD, 4);
        const double t4 = pp + 1;
        const double t5 = pp - 1;
        const double t6 = t4*t5;
        const double t7 = std::pow(ss, -5);
        const double t8 = 15*t7;
        const double t9 = std::pow(pp, 4);
        const double t10 = 35*t9;
        const double t11 = -35*t0 + t10;
        const double t12 = 3*t7;
        const double t13 = std::pow(DD, -3);
        const double t14 = t13*t6;
        const double t15 = std::pow(DD, -5);
        const double t16 = std::pow(ss, -4);
        const double t17 = (3.0/2.0)*t16;
        const double t18 = t1 - 4;
        const double t19 = std::pow(t4, 2)*std::pow(t5, 2);
        const double t20 = pp/std::pow(DD, 2);
        const double t21 = -45*t0 + t10 + 12;
        const double t22 = std::pow(ss, -3);
        const double t23 = pp/std::pow(DD, 6);
        const double t24 = -40*t0 + t10 + 8;
        const double t25 = t19/DD;
        const double t26 = -25*t0 + t10 + 2;
        const double t27 = (3.0/4.0)*t22;
        const double t28 = std::pow(DD, -7);
        const double t29 = (15.0/8.0)/std::pow(ss, 2);
        const double t30 = 120*t0;
        const double t31 = std::pow(pp, 6);
        const double t32 = std::pow(pp, 8);
        const double t33 = (105.0/16.0)/ss;
        o_mu_mu_mu_mu_mu[i] = t2*t3*t6*t8;
        o_mu_mu_mu_mu_sigma[i] = t12*t14*(t11 + 4);
        o_mu_mu_mu_mu_nu[i] = -t15*t17*(-30*t0 + t10 + 3);
        o_mu_mu_mu_sigma_sigma[i] = t18*t19*t20*t8;
        o_mu_mu_mu_sigma_nu[i] = -t17*t21*t3;
        o_mu_mu_mu_nu_nu[i] = (15.0/4.0)*t2*t22*t23;
        o_mu_mu_sigma_sigma_sigma[i] = t12*t24*t25;
        o_mu_mu_sigma_sigma_nu[i] = -t14*t17*t26;
        o_mu_mu_sigma_nu_nu[i] = t15*t26*t27;
        o_mu_mu_nu_nu_nu[i] = -t28*t29*(t1 - 1);
        o_mu_sigma_sigma_sigma_sigma[i] = pp*t12*(-t30 - 130*t31 + 35*t32 + 183*t9 + 40);
        o_mu_sigma_sigma_sigma_nu[i] = -t17*t20*t24*t6;
        o_mu_sigma_sigma_nu_nu[i] = t27*t3*(t11 + 6);
        o_mu_sigma_nu_nu_nu[i] = -t23*t29*(t1 - 2);
        o_mu_nu_nu_nu_nu[i] = pp*t33/std::pow(DD, 8);
        o_sigma_sigma_sigma_sigma_sigma[i] = t12*(35*DD*std::pow(pp, 10) + DD*t30 + 255*DD*t31 - 150*DD*t32 - 220*DD*t9 - 8);
        o_sigma_sigma_sigma_sigma_nu[i] = -15.0/2.0*t0*t16*t18*t25;
        o_sigma_sigma_sigma_nu_nu[i] = t0*t13*t21*t27;
        o_sigma_sigma_nu_nu_nu[i] = -t0*t15*t2*t29;
        o_sigma_nu_nu_nu_nu[i] = t0*t28*t33;
        // -(32 B5 D^9 + 105)/(32 D^9) split, so that D^9 cannot overflow
        // far in the tail (|z| = 1e300) and give inf/inf
        o_nu_nu_nu_nu_nu[i] = -B5 - (105.0/32.0)*t28/std::pow(DD, 2);
  }
  return List::create(Named("mu_mu_mu_mu_mu") = o_mu_mu_mu_mu_mu,
                      Named("mu_mu_mu_mu_sigma") = o_mu_mu_mu_mu_sigma,
                      Named("mu_mu_mu_mu_nu") = o_mu_mu_mu_mu_nu,
                      Named("mu_mu_mu_sigma_sigma") = o_mu_mu_mu_sigma_sigma,
                      Named("mu_mu_mu_sigma_nu") = o_mu_mu_mu_sigma_nu,
                      Named("mu_mu_mu_nu_nu") = o_mu_mu_mu_nu_nu,
                      Named("mu_mu_sigma_sigma_sigma") = o_mu_mu_sigma_sigma_sigma,
                      Named("mu_mu_sigma_sigma_nu") = o_mu_mu_sigma_sigma_nu,
                      Named("mu_mu_sigma_nu_nu") = o_mu_mu_sigma_nu_nu,
                      Named("mu_mu_nu_nu_nu") = o_mu_mu_nu_nu_nu,
                      Named("mu_sigma_sigma_sigma_sigma") = o_mu_sigma_sigma_sigma_sigma,
                      Named("mu_sigma_sigma_sigma_nu") = o_mu_sigma_sigma_sigma_nu,
                      Named("mu_sigma_sigma_nu_nu") = o_mu_sigma_sigma_nu_nu,
                      Named("mu_sigma_nu_nu_nu") = o_mu_sigma_nu_nu_nu,
                      Named("mu_nu_nu_nu_nu") = o_mu_nu_nu_nu_nu,
                      Named("sigma_sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma_sigma,
                      Named("sigma_sigma_sigma_sigma_nu") = o_sigma_sigma_sigma_sigma_nu,
                      Named("sigma_sigma_sigma_nu_nu") = o_sigma_sigma_sigma_nu_nu,
                      Named("sigma_sigma_nu_nu_nu") = o_sigma_sigma_nu_nu_nu,
                      Named("sigma_nu_nu_nu_nu") = o_sigma_nu_nu_nu_nu,
                      Named("nu_nu_nu_nu_nu") = o_nu_nu_nu_nu_nu);
}
