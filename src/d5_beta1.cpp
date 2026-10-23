#include <Rcpp.h>
#include <cmath>
#include "d7_par.h"
using namespace Rcpp;

// Fifth-order derivatives of the beta1 log-density in the parameters
// (mu, phi), closed form. Derived in sympy from the log-density and
// common-subexpression eliminated; written out per component. Same
// argument recycling as the deriv4 kernel: a scalar parameter is read
// once, a vector is read at i. Components are named as the deriv4 kernel.

// [[Rcpp::export]]
List beta1_deriv5_cpp(NumericVector y, NumericVector mu, NumericVector phi, int threads = 1) {
    int n = y.size();
    NumericVector o0(n), o1(n), o2(n), o3(n), o4(n), o5(n);
    bool mu_s = (mu.size() == 1);
    bool phi_s = (phi.size() == 1);
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        // BEGIN-BODY
        double m = mu_s ? mu[0] : mu[i];
        double p = phi_s ? phi[0] : phi[i];
        double t0 = std::pow(p, 5);
        double t1 = m*p + 1;
        double t2 = R::psigamma(t1, 4);
        double t3 = m - 1;
        double t4 = -p*t3 + 1;
        double t5 = R::psigamma(t4, 4);
        double t6 = std::pow(m, 5);
        double t7 = std::pow(t3, 5);
        double t8 = std::pow(m, 4);
        double t9 = std::pow(m, 2);
        double t10 = std::pow(m, 3);
        double t11 = R::psigamma(t1, 3);
        double t12 = std::pow(p, 3);
        double t13 = 4*t12;
        double t14 = std::pow(p, 4);
        double t15 = R::psigamma(t4, 3);
        double t16 = R::psigamma(t1, 2);
        double t17 = 6*p;
        double t18 = std::pow(p, 2);
        double t19 = 6*m;
        double t20 = R::psigamma(t4, 2);
        double t21 = std::pow(t3, 2);
        double t22 = std::pow(t3, 3);
        o0[i] = -t0*t2 + t0*t5 + 24*(2*m - 1)*(-3*m - 2*t10 + t8 + 4*t9 + 1)/(t6*t7);
        o1[i] = -m*t14*t2 - t11*t13 - t13*t15 + t14*t3*t5;
        o2[i] = 6*p*t20 - t11*t18*t19 - t12*t2*t9 + t12*t21*t5 - 6*t15*t18*t3 - t16*t17;
        o3[i] = -t10*t18*t2 - t11*t17*t9 - t15*t17*t21 - t16*t19 + t18*t22*t5 + 6*t20*t3;
        o4[i] = -p*t2*t8 + p*std::pow(t3, 4)*t5 - 4*t10*t11 - 4*t15*t22;
        o5[i] = -t2*t6 + t5*t7 + R::psigamma(p, 4) + 48/t0;
        // END-BODY
    });
    return List::create(Named("mu_mu_mu_mu_mu") = o0,
                        Named("mu_mu_mu_mu_phi") = o1,
                        Named("mu_mu_mu_phi_phi") = o2,
                        Named("mu_mu_phi_phi_phi") = o3,
                        Named("mu_phi_phi_phi_phi") = o4,
                        Named("phi_phi_phi_phi_phi") = o5);
}

