#include <Rcpp.h>
#include <cmath>
#include "d7_par.h"
using namespace Rcpp;

// Fifth-order derivatives of the invgauss1 log-density in the parameters
// (mu, phi), closed form. Derived in sympy from the log-density and
// common-subexpression eliminated; written out per component. Same
// argument recycling as the deriv4 kernel: a scalar parameter is read
// once, a vector is read at i. Components are named as the deriv4 kernel.

// [[Rcpp::export]]
List invgauss1_deriv5_cpp(NumericVector y, NumericVector mu, NumericVector phi, int threads = 1) {
    int n = y.size();
    NumericVector o0(n), o1(n), o2(n), o3(n), o4(n), o5(n);
    bool mu_s = (mu.size() == 1);
    bool phi_s = (phi.size() == 1);
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        // BEGIN-BODY
        double m = mu_s ? mu[0] : mu[i];
        double p = phi_s ? phi[0] : phi[i];
        double yi = y[i];
        double t0 = -3*yi;
        double t1 = 2*m;
        double t2 = std::pow(m, 2);
        o0[i] = -120*(m + t0)/(std::pow(m, 7)*p);
        o1[i] = -12*(t1 - 5*yi)/(std::pow(m, 6)*std::pow(p, 2));
        o2[i] = -12*(m - 2*yi)/(std::pow(m, 5)*std::pow(p, 3));
        o3[i] = -6*(t0 + t1)/(std::pow(m, 4)*std::pow(p, 4));
        o4[i] = -24*(m - yi)/(std::pow(m, 3)*std::pow(p, 5));
        o5[i] = -12*(10*m*yi + p*t2*yi - 5*t2 - 5*std::pow(yi, 2))/(std::pow(p, 6)*t2*yi);
        // END-BODY
    });
    return List::create(Named("mu_mu_mu_mu_mu") = o0,
                        Named("mu_mu_mu_mu_phi") = o1,
                        Named("mu_mu_mu_phi_phi") = o2,
                        Named("mu_mu_phi_phi_phi") = o3,
                        Named("mu_phi_phi_phi_phi") = o4,
                        Named("phi_phi_phi_phi_phi") = o5);
}

