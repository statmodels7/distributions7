#include <Rcpp.h>
#include <cmath>
#include "d7_par.h"
using namespace Rcpp;

// Fifth-order derivatives of the lognormal1 log-density in the parameters
// (mu, sigma2), closed form. Derived in sympy from the log-density and
// common-subexpression eliminated; written out per component. Same
// argument recycling as the deriv4 kernel: a scalar parameter is read
// once, a vector is read at i. Components are named as the deriv4 kernel.

// [[Rcpp::export]]
List lognormal1_deriv5_cpp(NumericVector y, NumericVector mu, NumericVector sigma2, int threads = 1) {
    int n = y.size();
    NumericVector o0(n), o1(n), o2(n), o3(n), o4(n), o5(n);
    bool mu_s = (mu.size() == 1);
    bool sigma2_s = (sigma2.size() == 1);
    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        // BEGIN-BODY
        double m = mu_s ? mu[0] : mu[i];
        double s2 = sigma2_s ? sigma2[0] : sigma2[i];
        double yi = y[i];
        double t0 = std::pow(s2, -5);
        double t1 = std::log(yi);
        double t2 = std::pow(s2, -6);
        o0[i] = 0;
        o1[i] = 0;
        o2[i] = 0;
        o3[i] = 6/std::pow(s2, 4);
        o4[i] = -24*m*t0 + 24*t0*t1;
        o5[i] = -120*m*t1*t2 + 60*std::pow(t1, 2)*t2 - 12*t2*(-5*std::pow(m, 2) + s2);
        // END-BODY
    });
    return List::create(Named("mu_mu_mu_mu_mu") = o0,
                        Named("mu_mu_mu_mu_sigma2") = o1,
                        Named("mu_mu_mu_sigma2_sigma2") = o2,
                        Named("mu_mu_sigma2_sigma2_sigma2") = o3,
                        Named("mu_sigma2_sigma2_sigma2_sigma2") = o4,
                        Named("sigma2_sigma2_sigma2_sigma2_sigma2") = o5);
}

