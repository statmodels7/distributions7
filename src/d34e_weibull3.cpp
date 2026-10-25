// The third and fourth derivatives of the expected information of the
// weibull3, on the parameter scale, in closed form. Each order is its own
// kernel. The closed form of the expected information is differentiated
// symbolically, and every polygamma is written either as its monomial minus
// a remainder R_k of psi_diff.h (arguments of at least 1) or as a regular
// part psi^(k)(1 + x) plus its pole (arguments below 1), so the monomial
// parts cancel exactly and no large term is ever summed with its opposite.

#include <Rcpp.h>
#include "d7_par.h"
#include "psi_diff.h"
#include "d34e_run.h"
using namespace Rcpp;

namespace {

// r_4(a) = -6/a^4 - psi^(4)(a)   (the form of d34e_gamma2.cpp)
inline double d34_r4(double a) {
    if (a >= 50.0) {
        const double u = 1.0 / a, u2 = u * u;
        return u * u * u * u * u *
               (12.0 + u * (10.0 + u2 * (-7.0 + u2 * (12.0 + u2 * (-33.0)))));
    }
    return -6.0 / (a * a * a * a) - R::psigamma(a, 4);
}

// r_5(a) = 24/a^5 - psi^(5)(a)
inline double d34_r5(double a) {
    if (a >= 50.0) {
        const double u = 1.0 / a, u2 = u * u;
        return -u * u * u * u * u * u *
               (60.0 + u * (60.0 + u2 * (-56.0 + u2 * (120.0 + u2 * (-396.0 +
                   u2 * 1820.0)))));
    }
    return 24.0 / (a * a * a * a * a) - R::psigamma(a, 5);
}

}  // namespace

// [[Rcpp::export]]
List weibull3_dexpected3_cpp(NumericVector y, NumericVector mean, NumericVector sigma, int threads = 1) {
    d34e::Par A0(mean), A1(sigma);
    return d34e::run(y.size(), threads, d7::kMinCostly, 12,
                     [&](std::size_t i, d34e::Row o) {
        const double m = A0[i], s = A1[i];
        const double X0 = 1 + 1.0/s;
        const double F0_0 = R::digamma(X0);
        const double F0_1 = R::psigamma(X0, 1);
        const double F0_2 = R::psigamma(X0, 2);
        const double F0_3 = R::psigamma(X0, 3);
        o[0] = 24*pow(s, 2)/pow(m, 5);  // mean_mean_mean_mean_mean
        o[1] = -12*s/pow(m, 4);  // mean_mean_mean_mean_sigma
        o[2] = 4/pow(m, 3);  // mean_mean_mean_sigma_sigma
        o[3] = 0;  // mean_mean_sigma_sigma_sigma
        o[4] = 0;  // sigma_sigma_mean_mean_mean
        o[5] = 0;  // sigma_sigma_mean_mean_sigma
        o[6] = 0;  // sigma_sigma_mean_sigma_sigma
        o[7] = -F0_1*(-8*F0_1/pow(s, 5) + 6*(0.84556867019693428 - 2*F0_0)/pow(s, 4))/pow(s, 2) - F0_1*(-12*F0_1/pow(s, 5) - 2*F0_2/pow(s, 6) + 6*(0.84556867019693428 - 2*F0_0)/pow(s, 4))/pow(s, 2) + 2*F0_1*(2*F0_1/pow(s, 4) - 2*(0.84556867019693428 - 2*F0_0)/pow(s, 3))/pow(s, 3) - 20*F0_1*(0.84556867019693428 - 2*F0_0)/pow(s, 6) - F0_2*(-2*F0_1/pow(s, 6) - (2*F0_1/pow(s, 4) - 2*(0.84556867019693428 - 2*F0_0)/pow(s, 3))/pow(s, 2) + 4*(0.84556867019693428 - 2*F0_0)/pow(s, 5))/pow(s, 2) - 6*F0_2*(0.84556867019693428 - 2*F0_0)/pow(s, 7) - F0_3*(0.84556867019693428 - 2*F0_0)/pow(s, 8) - 24*(-pow(F0_0 - 0.42278433509846714, 2) - 1.6449340668482264)/pow(s, 5);  // sigma_sigma_sigma_sigma_sigma
        o[8] = -6*(0.42278433509846714 - F0_0)/pow(m, 4);  // mean_sigma_mean_mean_mean
        o[9] = 2*F0_1/(pow(m, 3)*pow(s, 2));  // mean_sigma_mean_mean_sigma
        o[10] = 2*F0_1/(pow(m, 2)*pow(s, 3)) + F0_2/(pow(m, 2)*pow(s, 4));  // mean_sigma_mean_sigma_sigma
        o[11] = 6*F0_1/(m*pow(s, 4)) + 6*F0_2/(m*pow(s, 5)) + F0_3/(m*pow(s, 6));  // mean_sigma_sigma_sigma_sigma
    });
}

// [[Rcpp::export]]
List weibull3_dexpected4_cpp(NumericVector y, NumericVector mean, NumericVector sigma, int threads = 1) {
    d34e::Par A0(mean), A1(sigma);
    return d34e::run(y.size(), threads, d7::kMinCostly, 15,
                     [&](std::size_t i, d34e::Row o) {
        const double m = A0[i], s = A1[i];
        const double X0 = 1 + 1.0/s;
        const double F0_0 = R::digamma(X0);
        const double F0_1 = R::psigamma(X0, 1);
        const double F0_2 = R::psigamma(X0, 2);
        const double F0_3 = R::psigamma(X0, 3);
        const double F0_4 = R::psigamma(X0, 4);
        o[0] = -120*pow(s, 2)/pow(m, 6);  // mean_mean_mean_mean_mean_mean
        o[1] = 48*s/pow(m, 5);  // mean_mean_mean_mean_mean_sigma
        o[2] = -12/pow(m, 4);  // mean_mean_mean_mean_sigma_sigma
        o[3] = 0;  // mean_mean_mean_sigma_sigma_sigma
        o[4] = 0;  // mean_mean_sigma_sigma_sigma_sigma
        o[5] = 0;  // sigma_sigma_mean_mean_mean_mean
        o[6] = 0;  // sigma_sigma_mean_mean_mean_sigma
        o[7] = 0;  // sigma_sigma_mean_mean_sigma_sigma
        o[8] = 0;  // sigma_sigma_mean_sigma_sigma_sigma
        o[9] = -F0_1*(40*F0_1/pow(s, 6) - 24*(0.84556867019693428 - 2*F0_0)/pow(s, 5))/pow(s, 2) - F0_1*(60*F0_1/pow(s, 6) + 12*F0_2/pow(s, 7) - 24*(0.84556867019693428 - 2*F0_0)/pow(s, 5))/pow(s, 2) - F0_1*(72*F0_1/pow(s, 6) + 24*F0_2/pow(s, 7) + 2*F0_3/pow(s, 8) - 24*(0.84556867019693428 - 2*F0_0)/pow(s, 5))/pow(s, 2) + 4*F0_1*(-8*F0_1/pow(s, 5) + 6*(0.84556867019693428 - 2*F0_0)/pow(s, 4))/pow(s, 3) + 2*F0_1*(-12*F0_1/pow(s, 5) - 2*F0_2/pow(s, 6) + 6*(0.84556867019693428 - 2*F0_0)/pow(s, 4))/pow(s, 3) - 6*F0_1*(2*F0_1/pow(s, 4) - 2*(0.84556867019693428 - 2*F0_0)/pow(s, 3))/pow(s, 4) + 120*F0_1*(0.84556867019693428 - 2*F0_0)/pow(s, 7) - F0_2*(12*F0_1/pow(s, 7) - (-8*F0_1/pow(s, 5) + 6*(0.84556867019693428 - 2*F0_0)/pow(s, 4))/pow(s, 2) + 2*(2*F0_1/pow(s, 4) - 2*(0.84556867019693428 - 2*F0_0)/pow(s, 3))/pow(s, 3) - 20*(0.84556867019693428 - 2*F0_0)/pow(s, 6))/pow(s, 2) - F0_2*(24*F0_1/pow(s, 7) + 4*F0_2/pow(s, 8) - (-8*F0_1/pow(s, 5) + 6*(0.84556867019693428 - 2*F0_0)/pow(s, 4))/pow(s, 2) - (-12*F0_1/pow(s, 5) - 2*F0_2/pow(s, 6) + 6*(0.84556867019693428 - 2*F0_0)/pow(s, 4))/pow(s, 2) + 2*(2*F0_1/pow(s, 4) - 2*(0.84556867019693428 - 2*F0_0)/pow(s, 3))/pow(s, 3) - 20*(0.84556867019693428 - 2*F0_0)/pow(s, 6))/pow(s, 2) + 2*F0_2*(-2*F0_1/pow(s, 6) - (2*F0_1/pow(s, 4) - 2*(0.84556867019693428 - 2*F0_0)/pow(s, 3))/pow(s, 2) + 4*(0.84556867019693428 - 2*F0_0)/pow(s, 5))/pow(s, 3) + 42*F0_2*(0.84556867019693428 - 2*F0_0)/pow(s, 8) - F0_3*(2*F0_1/pow(s, 8) - (-2*F0_1/pow(s, 6) - (2*F0_1/pow(s, 4) - 2*(0.84556867019693428 - 2*F0_0)/pow(s, 3))/pow(s, 2) + 4*(0.84556867019693428 - 2*F0_0)/pow(s, 5))/pow(s, 2) - 6*(0.84556867019693428 - 2*F0_0)/pow(s, 7))/pow(s, 2) + 8*F0_3*(0.84556867019693428 - 2*F0_0)/pow(s, 9) + F0_4*(0.84556867019693428 - 2*F0_0)/pow(s, 10) + 120*(-pow(F0_0 - 0.42278433509846714, 2) - 1.6449340668482264)/pow(s, 6);  // sigma_sigma_sigma_sigma_sigma_sigma
        o[10] = 24*(0.42278433509846714 - F0_0)/pow(m, 5);  // mean_sigma_mean_mean_mean_mean
        o[11] = -6*F0_1/(pow(m, 4)*pow(s, 2));  // mean_sigma_mean_mean_mean_sigma
        o[12] = -4*F0_1/(pow(m, 3)*pow(s, 3)) - 2*F0_2/(pow(m, 3)*pow(s, 4));  // mean_sigma_mean_mean_sigma_sigma
        o[13] = -6*F0_1/(pow(m, 2)*pow(s, 4)) - 6*F0_2/(pow(m, 2)*pow(s, 5)) - F0_3/(pow(m, 2)*pow(s, 6));  // mean_sigma_mean_sigma_sigma_sigma
        o[14] = -24*F0_1/(m*pow(s, 5)) - 36*F0_2/(m*pow(s, 6)) - 12*F0_3/(m*pow(s, 7)) - F0_4/(m*pow(s, 8));  // mean_sigma_sigma_sigma_sigma_sigma
    });
}
