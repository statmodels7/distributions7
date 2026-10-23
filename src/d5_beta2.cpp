#include <Rcpp.h>
#include <cmath>
#include "d7_par.h"
using namespace Rcpp;

// Fifth-order derivatives of the beta2 log-density in the parameters
// (alpha, beta), closed form. Derived in sympy from the log-density and
// common-subexpression eliminated; written out per component. Same
// argument recycling as the deriv4 kernel: a scalar parameter is read
// once, a vector is read at i. Components are named as the deriv4 kernel.

// [[Rcpp::export]]
List beta2_deriv5_cpp(NumericVector y, NumericVector alpha, NumericVector beta) {
    int n = y.size();
    NumericVector o0(n), o1(n), o2(n), o3(n), o4(n), o5(n);
    bool alpha_s = (alpha.size() == 1);
    bool beta_s = (beta.size() == 1);
    for (int i = 0; i < n; ++i) {
        // BEGIN-BODY
        double a = alpha_s ? alpha[0] : alpha[i];
        double b = beta_s ? beta[0] : beta[i];
        double t0 = R::psigamma(a + b, 4);
        double t1 = -t0;
        o0[i] = -t1 - R::psigamma(a, 4);
        o1[i] = t0;
        o2[i] = t0;
        o3[i] = t0;
        o4[i] = t0;
        o5[i] = -t1 - R::psigamma(b, 4);
        // END-BODY
    }
    return List::create(Named("alpha_alpha_alpha_alpha_alpha") = o0,
                        Named("alpha_alpha_alpha_alpha_beta") = o1,
                        Named("alpha_alpha_alpha_beta_beta") = o2,
                        Named("alpha_alpha_beta_beta_beta") = o3,
                        Named("alpha_beta_beta_beta_beta") = o4,
                        Named("beta_beta_beta_beta_beta") = o5);
}

