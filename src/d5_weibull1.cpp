#include <Rcpp.h>
#include <cmath>
#include "d7_par.h"
using namespace Rcpp;

// Fifth-order derivatives of the weibull1 log-density in the parameters
// (mu, sigma), closed form. Derived in sympy from the log-density and
// common-subexpression eliminated; written out per component. Same
// argument recycling as the deriv4 kernel: a scalar parameter is read
// once, a vector is read at i. Components are named as the deriv4 kernel.

// [[Rcpp::export]]
List weibull1_deriv5_cpp(NumericVector y, NumericVector mu, NumericVector sigma, int threads = 1) {
    int n = y.size();
    NumericVector o0(n), o1(n), o2(n), o3(n), o4(n), o5(n);
    bool mu_s = (mu.size() == 1);
    bool sigma_s = (sigma.size() == 1);
    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        // BEGIN-BODY
        double m = mu_s ? mu[0] : mu[i];
        double s = sigma_s ? sigma[0] : sigma[i];
        double yi = y[i];
        double t0 = s/std::pow(m, 5);
        double t1 = 1.0/m;
        double t2 = std::log(t1*yi);
        double t3 = s*t2;
        double t4 = std::exp(t3);
        double t5 = std::pow(m, -4);
        double t6 = std::pow(s, 3);
        double t7 = std::pow(s, 2);
        double t8 = 6*t2;
        double t9 = 6*t3 + 6;
        double t10 = std::pow(t2, 2);
        double t11 = s*t10;
        double t12 = t10*t7;
        o0[i] = t0*t4*(s + 1)*(s + 2)*(s + 3)*(s + 4) - 24*t0;
        o1[i] = -t4*t5*(std::pow(s, 4)*t2 + 22*s + 11*t2*t7 + t6*t8 + 4*t6 + 18*t7 + t9) + 6*t5;
        o2[i] = t4*(6*s + t10*t6 + 2*t11 + 3*t12 + 4*t2 + 12*t3 + t7*t8 + 6)/std::pow(m, 3);
        o3[i] = -t2*t4*(t11 + t12 + 3*t2 + t9)/std::pow(m, 2);
        o4[i] = t1*std::pow(t2, 3)*t4*(t3 + 4);
        o5[i] = -std::pow(t2, 5)*t4 + 24/std::pow(s, 5);
        // END-BODY
    });
    return List::create(Named("mu_mu_mu_mu_mu") = o0,
                        Named("mu_mu_mu_mu_sigma") = o1,
                        Named("mu_mu_mu_sigma_sigma") = o2,
                        Named("mu_mu_sigma_sigma_sigma") = o3,
                        Named("mu_sigma_sigma_sigma_sigma") = o4,
                        Named("sigma_sigma_sigma_sigma_sigma") = o5);
}

