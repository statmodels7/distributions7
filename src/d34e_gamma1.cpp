// The third and fourth derivatives of the expected information of the
// gamma1, on the parameter scale, in closed form. Each order is its own
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
List gamma1_dexpected3_cpp(NumericVector y, NumericVector mu, NumericVector phi, int threads = 1) {
    d34e::Par A0(mu), A1(phi);
    return d34e::run(y.size(), threads, d7::kMinCostly, 12,
                     [&](std::size_t i, d34e::Row o) {
        const double m = A0[i], p = A1[i];
        const double X0 = 1.0/p;
        switch ((X0 < 1.0 ? 1 : 0)) {
        case 0: {
            const double R0_1 = d7::psi1_rest(X0);
            const double R0_2 = d7::psi2_rest(X0);
            const double R0_3 = d7::psi3_rest(X0);
            const double R0_4 = d34_r4(X0);
            o[0] = (24/(pow(m, 5)*p));  // mu_mu_mu_mu_mu
            o[1] = (6/(pow(m, 4)*pow(p, 2)));  // mu_mu_mu_mu_phi
            o[2] = (4/(pow(m, 3)*pow(p, 3)));  // mu_mu_mu_phi_phi
            o[3] = (6/(pow(m, 2)*pow(p, 4)));  // mu_mu_phi_phi_phi
            o[4] = 0.0;  // phi_phi_mu_mu_mu
            o[5] = 0.0;  // phi_phi_mu_mu_phi
            o[6] = 0.0;  // phi_phi_mu_phi_phi
            o[7] = (-120/pow(p, 7))*R0_1 + (-90/pow(p, 8))*R0_2 + (-18/pow(p, 9))*R0_3 + (-1/pow(p, 10))*R0_4;  // phi_phi_phi_phi_phi
            o[8] = 0.0;  // mu_phi_mu_mu_mu
            o[9] = 0.0;  // mu_phi_mu_mu_phi
            o[10] = 0.0;  // mu_phi_mu_phi_phi
            o[11] = 0.0;  // mu_phi_phi_phi_phi
            break;
        }
        case 1: {
            const double Q0_1 = R::psigamma(1.0 + X0, 1);
            const double Q0_2 = R::psigamma(1.0 + X0, 2);
            const double Q0_3 = R::psigamma(1.0 + X0, 3);
            const double Q0_4 = R::psigamma(1.0 + X0, 4);
            o[0] = (24/(pow(m, 5)*p));  // mu_mu_mu_mu_mu
            o[1] = (6/(pow(m, 4)*pow(p, 2)));  // mu_mu_mu_mu_phi
            o[2] = (4/(pow(m, 3)*pow(p, 3)));  // mu_mu_mu_phi_phi
            o[3] = (6/(pow(m, 2)*pow(p, 4)));  // mu_mu_phi_phi_phi
            o[4] = 0.0;  // phi_phi_mu_mu_mu
            o[5] = 0.0;  // phi_phi_mu_mu_phi
            o[6] = 0.0;  // phi_phi_mu_phi_phi
            o[7] = (12*(2*p - 5)/pow(p, 6)) + (120/pow(p, 7))*Q0_1 + (90/pow(p, 8))*Q0_2 + (18/pow(p, 9))*Q0_3 + (pow(p, -10))*Q0_4;  // phi_phi_phi_phi_phi
            o[8] = 0.0;  // mu_phi_mu_mu_mu
            o[9] = 0.0;  // mu_phi_mu_mu_phi
            o[10] = 0.0;  // mu_phi_mu_phi_phi
            o[11] = 0.0;  // mu_phi_phi_phi_phi
            break;
        }
        }
    });
}

// [[Rcpp::export]]
List gamma1_dexpected4_cpp(NumericVector y, NumericVector mu, NumericVector phi, int threads = 1) {
    d34e::Par A0(mu), A1(phi);
    return d34e::run(y.size(), threads, d7::kMinCostly, 15,
                     [&](std::size_t i, d34e::Row o) {
        const double m = A0[i], p = A1[i];
        const double X0 = 1.0/p;
        switch ((X0 < 1.0 ? 1 : 0)) {
        case 0: {
            const double R0_1 = d7::psi1_rest(X0);
            const double R0_2 = d7::psi2_rest(X0);
            const double R0_3 = d7::psi3_rest(X0);
            const double R0_4 = d34_r4(X0);
            const double R0_5 = d34_r5(X0);
            o[0] = (-120/(pow(m, 6)*p));  // mu_mu_mu_mu_mu_mu
            o[1] = (-24/(pow(m, 5)*pow(p, 2)));  // mu_mu_mu_mu_mu_phi
            o[2] = (-12/(pow(m, 4)*pow(p, 3)));  // mu_mu_mu_mu_phi_phi
            o[3] = (-12/(pow(m, 3)*pow(p, 4)));  // mu_mu_mu_phi_phi_phi
            o[4] = (-24/(pow(m, 2)*pow(p, 5)));  // mu_mu_phi_phi_phi_phi
            o[5] = 0.0;  // phi_phi_mu_mu_mu_mu
            o[6] = 0.0;  // phi_phi_mu_mu_mu_phi
            o[7] = 0.0;  // phi_phi_mu_mu_phi_phi
            o[8] = 0.0;  // phi_phi_mu_phi_phi_phi
            o[9] = (840/pow(p, 8))*R0_1 + (840/pow(p, 9))*R0_2 + (252/pow(p, 10))*R0_3 + (28/pow(p, 11))*R0_4 + (pow(p, -12))*R0_5;  // phi_phi_phi_phi_phi_phi
            o[10] = 0.0;  // mu_phi_mu_mu_mu_mu
            o[11] = 0.0;  // mu_phi_mu_mu_mu_phi
            o[12] = 0.0;  // mu_phi_mu_mu_phi_phi
            o[13] = 0.0;  // mu_phi_mu_phi_phi_phi
            o[14] = 0.0;  // mu_phi_phi_phi_phi_phi
            break;
        }
        case 1: {
            const double Q0_1 = R::psigamma(1.0 + X0, 1);
            const double Q0_2 = R::psigamma(1.0 + X0, 2);
            const double Q0_3 = R::psigamma(1.0 + X0, 3);
            const double Q0_4 = R::psigamma(1.0 + X0, 4);
            const double Q0_5 = R::psigamma(1.0 + X0, 5);
            o[0] = (-120/(pow(m, 6)*p));  // mu_mu_mu_mu_mu_mu
            o[1] = (-24/(pow(m, 5)*pow(p, 2)));  // mu_mu_mu_mu_mu_phi
            o[2] = (-12/(pow(m, 4)*pow(p, 3)));  // mu_mu_mu_mu_phi_phi
            o[3] = (-12/(pow(m, 3)*pow(p, 4)));  // mu_mu_mu_phi_phi_phi
            o[4] = (-24/(pow(m, 2)*pow(p, 5)));  // mu_mu_phi_phi_phi_phi
            o[5] = 0.0;  // phi_phi_mu_mu_mu_mu
            o[6] = 0.0;  // phi_phi_mu_mu_mu_phi
            o[7] = 0.0;  // phi_phi_mu_mu_phi_phi
            o[8] = 0.0;  // phi_phi_mu_phi_phi_phi
            o[9] = (-120*(p - 3)/pow(p, 7)) + (-840/pow(p, 8))*Q0_1 + (-840/pow(p, 9))*Q0_2 + (-252/pow(p, 10))*Q0_3 + (-28/pow(p, 11))*Q0_4 + (-1/pow(p, 12))*Q0_5;  // phi_phi_phi_phi_phi_phi
            o[10] = 0.0;  // mu_phi_mu_mu_mu_mu
            o[11] = 0.0;  // mu_phi_mu_mu_mu_phi
            o[12] = 0.0;  // mu_phi_mu_mu_phi_phi
            o[13] = 0.0;  // mu_phi_mu_phi_phi_phi
            o[14] = 0.0;  // mu_phi_phi_phi_phi_phi
            break;
        }
        }
    });
}
