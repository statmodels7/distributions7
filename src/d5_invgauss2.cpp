#include <Rcpp.h>
#include <cmath>
#include "d7_par.h"
using namespace Rcpp;

// Fifth-order derivatives of the invgauss2 log-density in the parameters
// (mu, lambda), closed form. Derived in sympy from the log-density and
// common-subexpression eliminated; written out per component. Same
// argument recycling as the deriv4 kernel: a scalar parameter is read
// once, a vector is read at i. Components are named as the deriv4 kernel.

// [[Rcpp::export]]
List invgauss2_deriv5_cpp(NumericVector y, NumericVector mu, NumericVector lambda, int threads = 1) {
    int n = y.size();
    NumericVector o0(n), o1(n), o2(n), o3(n), o4(n), o5(n);
    bool mu_s = (mu.size() == 1);
    bool lambda_s = (lambda.size() == 1);
    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        // BEGIN-BODY
        double m = mu_s ? mu[0] : mu[i];
        double L = lambda_s ? lambda[0] : lambda[i];
        double yi = y[i];
        o0[i] = -120*L*(m - 3*yi)/std::pow(m, 7);
        o1[i] = 12*(2*m - 5*yi)/std::pow(m, 6);
        o2[i] = 0;
        o3[i] = 0;
        o4[i] = 0;
        o5[i] = 12/std::pow(L, 5);
        // END-BODY
    });
    return List::create(Named("mu_mu_mu_mu_mu") = o0,
                        Named("mu_mu_mu_mu_lambda") = o1,
                        Named("mu_mu_mu_lambda_lambda") = o2,
                        Named("mu_mu_lambda_lambda_lambda") = o3,
                        Named("mu_lambda_lambda_lambda_lambda") = o4,
                        Named("lambda_lambda_lambda_lambda_lambda") = o5);
}

