// The third and fourth derivatives of the expected information of the
// beta2, on the parameter scale, in closed form. Each order is its own
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
List beta2_dexpected3_cpp(NumericVector y, NumericVector alpha, NumericVector beta, int threads = 1) {
    d34e::Par A0(alpha), A1(beta);
    return d34e::run(y.size(), threads, d7::kMinCostly, 12,
                     [&](std::size_t i, d34e::Row o) {
        const double a = A0[i], b = A1[i];
        const double X0 = a;
        const double X1 = b;
        const double X2 = a + b;
        switch ((X0 < 1.0 ? 1 : 0) | (X1 < 1.0 ? 2 : 0) | (X2 < 1.0 ? 4 : 0)) {
        case 0: {
            const double R0_4 = d34_r4(X0);
            const double R1_4 = d34_r4(X1);
            const double R2_4 = d34_r4(X2);
            o[0] = (6*b*(2*a + b)*(2*pow(a, 2) + 2*a*b + pow(b, 2))/(pow(a, 4)*pow(a + b, 4))) + (1)*R0_4 + (-1)*R2_4;  // alpha_alpha_alpha_alpha_alpha
            o[1] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_alpha_alpha_alpha_beta
            o[2] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_alpha_alpha_beta_beta
            o[3] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_alpha_beta_beta_beta
            o[4] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // beta_beta_alpha_alpha_alpha
            o[5] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // beta_beta_alpha_alpha_beta
            o[6] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // beta_beta_alpha_beta_beta
            o[7] = (6*a*(a + 2*b)*(pow(a, 2) + 2*a*b + 2*pow(b, 2))/(pow(b, 4)*pow(a + b, 4))) + (1)*R1_4 + (-1)*R2_4;  // beta_beta_beta_beta_beta
            o[8] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_beta_alpha_alpha_alpha
            o[9] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_beta_alpha_alpha_beta
            o[10] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_beta_alpha_beta_beta
            o[11] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_beta_beta_beta_beta
            break;
        }
        case 4: {
            const double Q2_4 = R::psigamma(1.0 + X2, 4);
            const double R0_4 = d34_r4(X0);
            const double R1_4 = d34_r4(X1);
            o[0] = (6*(pow(a, 5) + 5*pow(a, 4)*b - 4*pow(a, 4) + 10*pow(a, 3)*pow(b, 2) + 10*pow(a, 2)*pow(b, 3) + 5*a*pow(b, 4) + pow(b, 5))/(pow(a, 4)*pow(a + b, 5))) + (1)*R0_4 + (1)*Q2_4;  // alpha_alpha_alpha_alpha_alpha
            o[1] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_alpha_alpha_alpha_beta
            o[2] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_alpha_alpha_beta_beta
            o[3] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_alpha_beta_beta_beta
            o[4] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // beta_beta_alpha_alpha_alpha
            o[5] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // beta_beta_alpha_alpha_beta
            o[6] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // beta_beta_alpha_beta_beta
            o[7] = (6*(pow(a, 5) + 5*pow(a, 4)*b + 10*pow(a, 3)*pow(b, 2) + 10*pow(a, 2)*pow(b, 3) + 5*a*pow(b, 4) + pow(b, 5) - 4*pow(b, 4))/(pow(b, 4)*pow(a + b, 5))) + (1)*R1_4 + (1)*Q2_4;  // beta_beta_beta_beta_beta
            o[8] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_beta_alpha_alpha_alpha
            o[9] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_beta_alpha_alpha_beta
            o[10] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_beta_alpha_beta_beta
            o[11] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_beta_beta_beta_beta
            break;
        }
        case 2: {
            const double Q1_4 = R::psigamma(1.0 + X1, 4);
            const double R0_4 = d34_r4(X0);
            const double R2_4 = d34_r4(X2);
            o[0] = (6*b*(2*a + b)*(2*pow(a, 2) + 2*a*b + pow(b, 2))/(pow(a, 4)*pow(a + b, 4))) + (1)*R0_4 + (-1)*R2_4;  // alpha_alpha_alpha_alpha_alpha
            o[1] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_alpha_alpha_alpha_beta
            o[2] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_alpha_alpha_beta_beta
            o[3] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_alpha_beta_beta_beta
            o[4] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // beta_beta_alpha_alpha_alpha
            o[5] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // beta_beta_alpha_alpha_beta
            o[6] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // beta_beta_alpha_beta_beta
            o[7] = (6*(4*pow(a, 4) + 16*pow(a, 3)*b + 24*pow(a, 2)*pow(b, 2) + 16*a*pow(b, 3) - pow(b, 5) + 4*pow(b, 4))/(pow(b, 5)*pow(a + b, 4))) + (-1)*Q1_4 + (-1)*R2_4;  // beta_beta_beta_beta_beta
            o[8] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_beta_alpha_alpha_alpha
            o[9] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_beta_alpha_alpha_beta
            o[10] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_beta_alpha_beta_beta
            o[11] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_beta_beta_beta_beta
            break;
        }
        case 6: {
            const double Q1_4 = R::psigamma(1.0 + X1, 4);
            const double Q2_4 = R::psigamma(1.0 + X2, 4);
            const double R0_4 = d34_r4(X0);
            o[0] = (6*(pow(a, 5) + 5*pow(a, 4)*b - 4*pow(a, 4) + 10*pow(a, 3)*pow(b, 2) + 10*pow(a, 2)*pow(b, 3) + 5*a*pow(b, 4) + pow(b, 5))/(pow(a, 4)*pow(a + b, 5))) + (1)*R0_4 + (1)*Q2_4;  // alpha_alpha_alpha_alpha_alpha
            o[1] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_alpha_alpha_alpha_beta
            o[2] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_alpha_alpha_beta_beta
            o[3] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_alpha_beta_beta_beta
            o[4] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // beta_beta_alpha_alpha_alpha
            o[5] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // beta_beta_alpha_alpha_beta
            o[6] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // beta_beta_alpha_beta_beta
            o[7] = (24*a*(pow(a, 4) + 5*pow(a, 3)*b + 10*pow(a, 2)*pow(b, 2) + 10*a*pow(b, 3) + 5*pow(b, 4))/(pow(b, 5)*pow(a + b, 5))) + (-1)*Q1_4 + (1)*Q2_4;  // beta_beta_beta_beta_beta
            o[8] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_beta_alpha_alpha_alpha
            o[9] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_beta_alpha_alpha_beta
            o[10] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_beta_alpha_beta_beta
            o[11] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_beta_beta_beta_beta
            break;
        }
        case 1: {
            const double Q0_4 = R::psigamma(1.0 + X0, 4);
            const double R1_4 = d34_r4(X1);
            const double R2_4 = d34_r4(X2);
            o[0] = (-6*(pow(a, 5) - 4*pow(a, 4) - 16*pow(a, 3)*b - 24*pow(a, 2)*pow(b, 2) - 16*a*pow(b, 3) - 4*pow(b, 4))/(pow(a, 5)*pow(a + b, 4))) + (-1)*Q0_4 + (-1)*R2_4;  // alpha_alpha_alpha_alpha_alpha
            o[1] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_alpha_alpha_alpha_beta
            o[2] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_alpha_alpha_beta_beta
            o[3] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_alpha_beta_beta_beta
            o[4] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // beta_beta_alpha_alpha_alpha
            o[5] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // beta_beta_alpha_alpha_beta
            o[6] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // beta_beta_alpha_beta_beta
            o[7] = (6*a*(a + 2*b)*(pow(a, 2) + 2*a*b + 2*pow(b, 2))/(pow(b, 4)*pow(a + b, 4))) + (1)*R1_4 + (-1)*R2_4;  // beta_beta_beta_beta_beta
            o[8] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_beta_alpha_alpha_alpha
            o[9] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_beta_alpha_alpha_beta
            o[10] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_beta_alpha_beta_beta
            o[11] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_beta_beta_beta_beta
            break;
        }
        case 5: {
            const double Q0_4 = R::psigamma(1.0 + X0, 4);
            const double Q2_4 = R::psigamma(1.0 + X2, 4);
            const double R1_4 = d34_r4(X1);
            o[0] = (24*b*(5*pow(a, 4) + 10*pow(a, 3)*b + 10*pow(a, 2)*pow(b, 2) + 5*a*pow(b, 3) + pow(b, 4))/(pow(a, 5)*pow(a + b, 5))) + (-1)*Q0_4 + (1)*Q2_4;  // alpha_alpha_alpha_alpha_alpha
            o[1] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_alpha_alpha_alpha_beta
            o[2] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_alpha_alpha_beta_beta
            o[3] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_alpha_beta_beta_beta
            o[4] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // beta_beta_alpha_alpha_alpha
            o[5] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // beta_beta_alpha_alpha_beta
            o[6] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // beta_beta_alpha_beta_beta
            o[7] = (6*(pow(a, 5) + 5*pow(a, 4)*b + 10*pow(a, 3)*pow(b, 2) + 10*pow(a, 2)*pow(b, 3) + 5*a*pow(b, 4) + pow(b, 5) - 4*pow(b, 4))/(pow(b, 4)*pow(a + b, 5))) + (1)*R1_4 + (1)*Q2_4;  // beta_beta_beta_beta_beta
            o[8] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_beta_alpha_alpha_alpha
            o[9] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_beta_alpha_alpha_beta
            o[10] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_beta_alpha_beta_beta
            o[11] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_beta_beta_beta_beta
            break;
        }
        case 3: {
            const double Q0_4 = R::psigamma(1.0 + X0, 4);
            const double Q1_4 = R::psigamma(1.0 + X1, 4);
            const double R2_4 = d34_r4(X2);
            o[0] = (-6*(pow(a, 5) - 4*pow(a, 4) - 16*pow(a, 3)*b - 24*pow(a, 2)*pow(b, 2) - 16*a*pow(b, 3) - 4*pow(b, 4))/(pow(a, 5)*pow(a + b, 4))) + (-1)*Q0_4 + (-1)*R2_4;  // alpha_alpha_alpha_alpha_alpha
            o[1] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_alpha_alpha_alpha_beta
            o[2] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_alpha_alpha_beta_beta
            o[3] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_alpha_beta_beta_beta
            o[4] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // beta_beta_alpha_alpha_alpha
            o[5] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // beta_beta_alpha_alpha_beta
            o[6] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // beta_beta_alpha_beta_beta
            o[7] = (6*(4*pow(a, 4) + 16*pow(a, 3)*b + 24*pow(a, 2)*pow(b, 2) + 16*a*pow(b, 3) - pow(b, 5) + 4*pow(b, 4))/(pow(b, 5)*pow(a + b, 4))) + (-1)*Q1_4 + (-1)*R2_4;  // beta_beta_beta_beta_beta
            o[8] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_beta_alpha_alpha_alpha
            o[9] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_beta_alpha_alpha_beta
            o[10] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_beta_alpha_beta_beta
            o[11] = (-6/pow(a + b, 4)) + (-1)*R2_4;  // alpha_beta_beta_beta_beta
            break;
        }
        case 7: {
            const double Q0_4 = R::psigamma(1.0 + X0, 4);
            const double Q1_4 = R::psigamma(1.0 + X1, 4);
            const double Q2_4 = R::psigamma(1.0 + X2, 4);
            o[0] = (24*b*(5*pow(a, 4) + 10*pow(a, 3)*b + 10*pow(a, 2)*pow(b, 2) + 5*a*pow(b, 3) + pow(b, 4))/(pow(a, 5)*pow(a + b, 5))) + (-1)*Q0_4 + (1)*Q2_4;  // alpha_alpha_alpha_alpha_alpha
            o[1] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_alpha_alpha_alpha_beta
            o[2] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_alpha_alpha_beta_beta
            o[3] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_alpha_beta_beta_beta
            o[4] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // beta_beta_alpha_alpha_alpha
            o[5] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // beta_beta_alpha_alpha_beta
            o[6] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // beta_beta_alpha_beta_beta
            o[7] = (24*a*(pow(a, 4) + 5*pow(a, 3)*b + 10*pow(a, 2)*pow(b, 2) + 10*a*pow(b, 3) + 5*pow(b, 4))/(pow(b, 5)*pow(a + b, 5))) + (-1)*Q1_4 + (1)*Q2_4;  // beta_beta_beta_beta_beta
            o[8] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_beta_alpha_alpha_alpha
            o[9] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_beta_alpha_alpha_beta
            o[10] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_beta_alpha_beta_beta
            o[11] = (-24/pow(a + b, 5)) + (1)*Q2_4;  // alpha_beta_beta_beta_beta
            break;
        }
        }
    });
}

// [[Rcpp::export]]
List beta2_dexpected4_cpp(NumericVector y, NumericVector alpha, NumericVector beta, int threads = 1) {
    d34e::Par A0(alpha), A1(beta);
    return d34e::run(y.size(), threads, d7::kMinCostly, 15,
                     [&](std::size_t i, d34e::Row o) {
        const double a = A0[i], b = A1[i];
        const double X0 = a;
        const double X1 = b;
        const double X2 = a + b;
        switch ((X0 < 1.0 ? 1 : 0) | (X1 < 1.0 ? 2 : 0) | (X2 < 1.0 ? 4 : 0)) {
        case 0: {
            const double R0_5 = d34_r5(X0);
            const double R1_5 = d34_r5(X1);
            const double R2_5 = d34_r5(X2);
            o[0] = (-24*b*(5*pow(a, 4) + 10*pow(a, 3)*b + 10*pow(a, 2)*pow(b, 2) + 5*a*pow(b, 3) + pow(b, 4))/(pow(a, 5)*pow(a + b, 5))) + (1)*R0_5 + (-1)*R2_5;  // alpha_alpha_alpha_alpha_alpha_alpha
            o[1] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_alpha_alpha_alpha_alpha_beta
            o[2] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_alpha_alpha_alpha_beta_beta
            o[3] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_alpha_alpha_beta_beta_beta
            o[4] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_alpha_beta_beta_beta_beta
            o[5] = (24/pow(a + b, 5)) + (-1)*R2_5;  // beta_beta_alpha_alpha_alpha_alpha
            o[6] = (24/pow(a + b, 5)) + (-1)*R2_5;  // beta_beta_alpha_alpha_alpha_beta
            o[7] = (24/pow(a + b, 5)) + (-1)*R2_5;  // beta_beta_alpha_alpha_beta_beta
            o[8] = (24/pow(a + b, 5)) + (-1)*R2_5;  // beta_beta_alpha_beta_beta_beta
            o[9] = (-24*a*(pow(a, 4) + 5*pow(a, 3)*b + 10*pow(a, 2)*pow(b, 2) + 10*a*pow(b, 3) + 5*pow(b, 4))/(pow(b, 5)*pow(a + b, 5))) + (1)*R1_5 + (-1)*R2_5;  // beta_beta_beta_beta_beta_beta
            o[10] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_beta_alpha_alpha_alpha_alpha
            o[11] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_beta_alpha_alpha_alpha_beta
            o[12] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_beta_alpha_alpha_beta_beta
            o[13] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_beta_alpha_beta_beta_beta
            o[14] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_beta_beta_beta_beta_beta
            break;
        }
        case 4: {
            const double Q2_5 = R::psigamma(1.0 + X2, 5);
            const double R0_5 = d34_r5(X0);
            const double R1_5 = d34_r5(X1);
            o[0] = (-24*(pow(a, 6) + 6*pow(a, 5)*b - 5*pow(a, 5) + 15*pow(a, 4)*pow(b, 2) + 20*pow(a, 3)*pow(b, 3) + 15*pow(a, 2)*pow(b, 4) + 6*a*pow(b, 5) + pow(b, 6))/(pow(a, 5)*pow(a + b, 6))) + (1)*R0_5 + (1)*Q2_5;  // alpha_alpha_alpha_alpha_alpha_alpha
            o[1] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_alpha_alpha_alpha_alpha_beta
            o[2] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_alpha_alpha_alpha_beta_beta
            o[3] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_alpha_alpha_beta_beta_beta
            o[4] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_alpha_beta_beta_beta_beta
            o[5] = (120/pow(a + b, 6)) + (1)*Q2_5;  // beta_beta_alpha_alpha_alpha_alpha
            o[6] = (120/pow(a + b, 6)) + (1)*Q2_5;  // beta_beta_alpha_alpha_alpha_beta
            o[7] = (120/pow(a + b, 6)) + (1)*Q2_5;  // beta_beta_alpha_alpha_beta_beta
            o[8] = (120/pow(a + b, 6)) + (1)*Q2_5;  // beta_beta_alpha_beta_beta_beta
            o[9] = (-24*(pow(a, 6) + 6*pow(a, 5)*b + 15*pow(a, 4)*pow(b, 2) + 20*pow(a, 3)*pow(b, 3) + 15*pow(a, 2)*pow(b, 4) + 6*a*pow(b, 5) + pow(b, 6) - 5*pow(b, 5))/(pow(b, 5)*pow(a + b, 6))) + (1)*R1_5 + (1)*Q2_5;  // beta_beta_beta_beta_beta_beta
            o[10] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_beta_alpha_alpha_alpha_alpha
            o[11] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_beta_alpha_alpha_alpha_beta
            o[12] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_beta_alpha_alpha_beta_beta
            o[13] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_beta_alpha_beta_beta_beta
            o[14] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_beta_beta_beta_beta_beta
            break;
        }
        case 2: {
            const double Q1_5 = R::psigamma(1.0 + X1, 5);
            const double R0_5 = d34_r5(X0);
            const double R2_5 = d34_r5(X2);
            o[0] = (-24*b*(5*pow(a, 4) + 10*pow(a, 3)*b + 10*pow(a, 2)*pow(b, 2) + 5*a*pow(b, 3) + pow(b, 4))/(pow(a, 5)*pow(a + b, 5))) + (1)*R0_5 + (-1)*R2_5;  // alpha_alpha_alpha_alpha_alpha_alpha
            o[1] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_alpha_alpha_alpha_alpha_beta
            o[2] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_alpha_alpha_alpha_beta_beta
            o[3] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_alpha_alpha_beta_beta_beta
            o[4] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_alpha_beta_beta_beta_beta
            o[5] = (24/pow(a + b, 5)) + (-1)*R2_5;  // beta_beta_alpha_alpha_alpha_alpha
            o[6] = (24/pow(a + b, 5)) + (-1)*R2_5;  // beta_beta_alpha_alpha_alpha_beta
            o[7] = (24/pow(a + b, 5)) + (-1)*R2_5;  // beta_beta_alpha_alpha_beta_beta
            o[8] = (24/pow(a + b, 5)) + (-1)*R2_5;  // beta_beta_alpha_beta_beta_beta
            o[9] = (-24*(5*pow(a, 5) + 25*pow(a, 4)*b + 50*pow(a, 3)*pow(b, 2) + 50*pow(a, 2)*pow(b, 3) + 25*a*pow(b, 4) - pow(b, 6) + 5*pow(b, 5))/(pow(b, 6)*pow(a + b, 5))) + (-1)*Q1_5 + (-1)*R2_5;  // beta_beta_beta_beta_beta_beta
            o[10] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_beta_alpha_alpha_alpha_alpha
            o[11] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_beta_alpha_alpha_alpha_beta
            o[12] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_beta_alpha_alpha_beta_beta
            o[13] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_beta_alpha_beta_beta_beta
            o[14] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_beta_beta_beta_beta_beta
            break;
        }
        case 6: {
            const double Q1_5 = R::psigamma(1.0 + X1, 5);
            const double Q2_5 = R::psigamma(1.0 + X2, 5);
            const double R0_5 = d34_r5(X0);
            o[0] = (-24*(pow(a, 6) + 6*pow(a, 5)*b - 5*pow(a, 5) + 15*pow(a, 4)*pow(b, 2) + 20*pow(a, 3)*pow(b, 3) + 15*pow(a, 2)*pow(b, 4) + 6*a*pow(b, 5) + pow(b, 6))/(pow(a, 5)*pow(a + b, 6))) + (1)*R0_5 + (1)*Q2_5;  // alpha_alpha_alpha_alpha_alpha_alpha
            o[1] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_alpha_alpha_alpha_alpha_beta
            o[2] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_alpha_alpha_alpha_beta_beta
            o[3] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_alpha_alpha_beta_beta_beta
            o[4] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_alpha_beta_beta_beta_beta
            o[5] = (120/pow(a + b, 6)) + (1)*Q2_5;  // beta_beta_alpha_alpha_alpha_alpha
            o[6] = (120/pow(a + b, 6)) + (1)*Q2_5;  // beta_beta_alpha_alpha_alpha_beta
            o[7] = (120/pow(a + b, 6)) + (1)*Q2_5;  // beta_beta_alpha_alpha_beta_beta
            o[8] = (120/pow(a + b, 6)) + (1)*Q2_5;  // beta_beta_alpha_beta_beta_beta
            o[9] = (-120*a*(a + 2*b)*(pow(a, 2) + a*b + pow(b, 2))*(pow(a, 2) + 3*a*b + 3*pow(b, 2))/(pow(b, 6)*pow(a + b, 6))) + (-1)*Q1_5 + (1)*Q2_5;  // beta_beta_beta_beta_beta_beta
            o[10] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_beta_alpha_alpha_alpha_alpha
            o[11] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_beta_alpha_alpha_alpha_beta
            o[12] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_beta_alpha_alpha_beta_beta
            o[13] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_beta_alpha_beta_beta_beta
            o[14] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_beta_beta_beta_beta_beta
            break;
        }
        case 1: {
            const double Q0_5 = R::psigamma(1.0 + X0, 5);
            const double R1_5 = d34_r5(X1);
            const double R2_5 = d34_r5(X2);
            o[0] = (24*(pow(a, 6) - 5*pow(a, 5) - 25*pow(a, 4)*b - 50*pow(a, 3)*pow(b, 2) - 50*pow(a, 2)*pow(b, 3) - 25*a*pow(b, 4) - 5*pow(b, 5))/(pow(a, 6)*pow(a + b, 5))) + (-1)*Q0_5 + (-1)*R2_5;  // alpha_alpha_alpha_alpha_alpha_alpha
            o[1] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_alpha_alpha_alpha_alpha_beta
            o[2] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_alpha_alpha_alpha_beta_beta
            o[3] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_alpha_alpha_beta_beta_beta
            o[4] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_alpha_beta_beta_beta_beta
            o[5] = (24/pow(a + b, 5)) + (-1)*R2_5;  // beta_beta_alpha_alpha_alpha_alpha
            o[6] = (24/pow(a + b, 5)) + (-1)*R2_5;  // beta_beta_alpha_alpha_alpha_beta
            o[7] = (24/pow(a + b, 5)) + (-1)*R2_5;  // beta_beta_alpha_alpha_beta_beta
            o[8] = (24/pow(a + b, 5)) + (-1)*R2_5;  // beta_beta_alpha_beta_beta_beta
            o[9] = (-24*a*(pow(a, 4) + 5*pow(a, 3)*b + 10*pow(a, 2)*pow(b, 2) + 10*a*pow(b, 3) + 5*pow(b, 4))/(pow(b, 5)*pow(a + b, 5))) + (1)*R1_5 + (-1)*R2_5;  // beta_beta_beta_beta_beta_beta
            o[10] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_beta_alpha_alpha_alpha_alpha
            o[11] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_beta_alpha_alpha_alpha_beta
            o[12] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_beta_alpha_alpha_beta_beta
            o[13] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_beta_alpha_beta_beta_beta
            o[14] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_beta_beta_beta_beta_beta
            break;
        }
        case 5: {
            const double Q0_5 = R::psigamma(1.0 + X0, 5);
            const double Q2_5 = R::psigamma(1.0 + X2, 5);
            const double R1_5 = d34_r5(X1);
            o[0] = (-120*b*(2*a + b)*(pow(a, 2) + a*b + pow(b, 2))*(3*pow(a, 2) + 3*a*b + pow(b, 2))/(pow(a, 6)*pow(a + b, 6))) + (-1)*Q0_5 + (1)*Q2_5;  // alpha_alpha_alpha_alpha_alpha_alpha
            o[1] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_alpha_alpha_alpha_alpha_beta
            o[2] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_alpha_alpha_alpha_beta_beta
            o[3] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_alpha_alpha_beta_beta_beta
            o[4] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_alpha_beta_beta_beta_beta
            o[5] = (120/pow(a + b, 6)) + (1)*Q2_5;  // beta_beta_alpha_alpha_alpha_alpha
            o[6] = (120/pow(a + b, 6)) + (1)*Q2_5;  // beta_beta_alpha_alpha_alpha_beta
            o[7] = (120/pow(a + b, 6)) + (1)*Q2_5;  // beta_beta_alpha_alpha_beta_beta
            o[8] = (120/pow(a + b, 6)) + (1)*Q2_5;  // beta_beta_alpha_beta_beta_beta
            o[9] = (-24*(pow(a, 6) + 6*pow(a, 5)*b + 15*pow(a, 4)*pow(b, 2) + 20*pow(a, 3)*pow(b, 3) + 15*pow(a, 2)*pow(b, 4) + 6*a*pow(b, 5) + pow(b, 6) - 5*pow(b, 5))/(pow(b, 5)*pow(a + b, 6))) + (1)*R1_5 + (1)*Q2_5;  // beta_beta_beta_beta_beta_beta
            o[10] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_beta_alpha_alpha_alpha_alpha
            o[11] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_beta_alpha_alpha_alpha_beta
            o[12] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_beta_alpha_alpha_beta_beta
            o[13] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_beta_alpha_beta_beta_beta
            o[14] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_beta_beta_beta_beta_beta
            break;
        }
        case 3: {
            const double Q0_5 = R::psigamma(1.0 + X0, 5);
            const double Q1_5 = R::psigamma(1.0 + X1, 5);
            const double R2_5 = d34_r5(X2);
            o[0] = (24*(pow(a, 6) - 5*pow(a, 5) - 25*pow(a, 4)*b - 50*pow(a, 3)*pow(b, 2) - 50*pow(a, 2)*pow(b, 3) - 25*a*pow(b, 4) - 5*pow(b, 5))/(pow(a, 6)*pow(a + b, 5))) + (-1)*Q0_5 + (-1)*R2_5;  // alpha_alpha_alpha_alpha_alpha_alpha
            o[1] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_alpha_alpha_alpha_alpha_beta
            o[2] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_alpha_alpha_alpha_beta_beta
            o[3] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_alpha_alpha_beta_beta_beta
            o[4] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_alpha_beta_beta_beta_beta
            o[5] = (24/pow(a + b, 5)) + (-1)*R2_5;  // beta_beta_alpha_alpha_alpha_alpha
            o[6] = (24/pow(a + b, 5)) + (-1)*R2_5;  // beta_beta_alpha_alpha_alpha_beta
            o[7] = (24/pow(a + b, 5)) + (-1)*R2_5;  // beta_beta_alpha_alpha_beta_beta
            o[8] = (24/pow(a + b, 5)) + (-1)*R2_5;  // beta_beta_alpha_beta_beta_beta
            o[9] = (-24*(5*pow(a, 5) + 25*pow(a, 4)*b + 50*pow(a, 3)*pow(b, 2) + 50*pow(a, 2)*pow(b, 3) + 25*a*pow(b, 4) - pow(b, 6) + 5*pow(b, 5))/(pow(b, 6)*pow(a + b, 5))) + (-1)*Q1_5 + (-1)*R2_5;  // beta_beta_beta_beta_beta_beta
            o[10] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_beta_alpha_alpha_alpha_alpha
            o[11] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_beta_alpha_alpha_alpha_beta
            o[12] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_beta_alpha_alpha_beta_beta
            o[13] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_beta_alpha_beta_beta_beta
            o[14] = (24/pow(a + b, 5)) + (-1)*R2_5;  // alpha_beta_beta_beta_beta_beta
            break;
        }
        case 7: {
            const double Q0_5 = R::psigamma(1.0 + X0, 5);
            const double Q1_5 = R::psigamma(1.0 + X1, 5);
            const double Q2_5 = R::psigamma(1.0 + X2, 5);
            o[0] = (-120*b*(2*a + b)*(pow(a, 2) + a*b + pow(b, 2))*(3*pow(a, 2) + 3*a*b + pow(b, 2))/(pow(a, 6)*pow(a + b, 6))) + (-1)*Q0_5 + (1)*Q2_5;  // alpha_alpha_alpha_alpha_alpha_alpha
            o[1] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_alpha_alpha_alpha_alpha_beta
            o[2] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_alpha_alpha_alpha_beta_beta
            o[3] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_alpha_alpha_beta_beta_beta
            o[4] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_alpha_beta_beta_beta_beta
            o[5] = (120/pow(a + b, 6)) + (1)*Q2_5;  // beta_beta_alpha_alpha_alpha_alpha
            o[6] = (120/pow(a + b, 6)) + (1)*Q2_5;  // beta_beta_alpha_alpha_alpha_beta
            o[7] = (120/pow(a + b, 6)) + (1)*Q2_5;  // beta_beta_alpha_alpha_beta_beta
            o[8] = (120/pow(a + b, 6)) + (1)*Q2_5;  // beta_beta_alpha_beta_beta_beta
            o[9] = (-120*a*(a + 2*b)*(pow(a, 2) + a*b + pow(b, 2))*(pow(a, 2) + 3*a*b + 3*pow(b, 2))/(pow(b, 6)*pow(a + b, 6))) + (-1)*Q1_5 + (1)*Q2_5;  // beta_beta_beta_beta_beta_beta
            o[10] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_beta_alpha_alpha_alpha_alpha
            o[11] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_beta_alpha_alpha_alpha_beta
            o[12] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_beta_alpha_alpha_beta_beta
            o[13] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_beta_alpha_beta_beta_beta
            o[14] = (120/pow(a + b, 6)) + (1)*Q2_5;  // alpha_beta_beta_beta_beta_beta
            break;
        }
        }
    });
}
