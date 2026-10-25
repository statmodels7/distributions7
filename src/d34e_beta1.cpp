// The third and fourth derivatives of the expected information of the
// beta1, on the parameter scale, in closed form. Each order is its own
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
List beta1_dexpected3_cpp(NumericVector y, NumericVector mu, NumericVector phi, int threads = 1) {
    d34e::Par A0(mu), A1(phi);
    return d34e::run(y.size(), threads, d7::kMinCostly, 12,
                     [&](std::size_t i, d34e::Row o) {
        const double m = A0[i], p = A1[i];
        const double X0 = m*p;
        const double X1 = p*(1 - m);
        const double X2 = p;
        switch ((X0 < 1.0 ? 1 : 0) | (X1 < 1.0 ? 2 : 0) | (X2 < 1.0 ? 4 : 0)) {
        case 0: {
            const double R0_2 = d7::psi2_rest(X0);
            const double R0_3 = d7::psi3_rest(X0);
            const double R0_4 = d34_r4(X0);
            const double R1_2 = d7::psi2_rest(X1);
            const double R1_3 = d7::psi3_rest(X1);
            const double R1_4 = d34_r4(X1);
            const double R2_4 = d34_r4(X2);
            o[0] = (-6*p*(2*m - 1)*(2*pow(m, 2) - 2*m + 1)/(pow(m, 4)*pow(m - 1, 4))) + (pow(p, 5))*R0_4 + (-pow(p, 5))*R1_4;  // mu_mu_mu_mu_mu
            o[1] = (2*(3*pow(m, 2) - 3*m + 1)/(pow(m, 3)*pow(m - 1, 3))) + (4*pow(p, 3))*R0_3 + (m*pow(p, 4))*R0_4 + (4*pow(p, 3))*R1_3 + (-pow(p, 4)*(m - 1))*R1_4;  // mu_mu_mu_mu_phi
            o[2] = (6*p)*R0_2 + (6*m*pow(p, 2))*R0_3 + (pow(m, 2)*pow(p, 3))*R0_4 + (-6*p)*R1_2 + (6*pow(p, 2)*(m - 1))*R1_3 + (-pow(p, 3)*pow(m - 1, 2))*R1_4;  // mu_mu_mu_phi_phi
            o[3] = (6*m)*R0_2 + (6*pow(m, 2)*p)*R0_3 + (pow(m, 3)*pow(p, 2))*R0_4 + (-6*(m - 1))*R1_2 + (6*p*pow(m - 1, 2))*R1_3 + (-pow(p, 2)*pow(m - 1, 3))*R1_4;  // mu_mu_phi_phi_phi
            o[4] = (6*p)*R0_2 + (6*m*pow(p, 2))*R0_3 + (pow(m, 2)*pow(p, 3))*R0_4 + (-6*p)*R1_2 + (6*pow(p, 2)*(m - 1))*R1_3 + (-pow(p, 3)*pow(m - 1, 2))*R1_4;  // phi_phi_mu_mu_mu
            o[5] = (6*m)*R0_2 + (6*pow(m, 2)*p)*R0_3 + (pow(m, 3)*pow(p, 2))*R0_4 + (-6*(m - 1))*R1_2 + (6*p*pow(m - 1, 2))*R1_3 + (-pow(p, 2)*pow(m - 1, 3))*R1_4;  // phi_phi_mu_mu_phi
            o[6] = (4*pow(m, 3))*R0_3 + (pow(m, 4)*p)*R0_4 + (4*pow(m - 1, 3))*R1_3 + (-p*pow(m - 1, 4))*R1_4;  // phi_phi_mu_phi_phi
            o[7] = (pow(m, 5))*R0_4 + (-pow(m - 1, 5))*R1_4 + (-1)*R2_4;  // phi_phi_phi_phi_phi
            o[8] = (3*pow(p, 3))*R0_3 + (m*pow(p, 4))*R0_4 + (3*pow(p, 3))*R1_3 + (-pow(p, 4)*(m - 1))*R1_4;  // mu_phi_mu_mu_mu
            o[9] = (4*p)*R0_2 + (5*m*pow(p, 2))*R0_3 + (pow(m, 2)*pow(p, 3))*R0_4 + (-4*p)*R1_2 + (5*pow(p, 2)*(m - 1))*R1_3 + (-pow(p, 3)*pow(m - 1, 2))*R1_4;  // mu_phi_mu_mu_phi
            o[10] = (4*m)*R0_2 + (5*pow(m, 2)*p)*R0_3 + (pow(m, 3)*pow(p, 2))*R0_4 + (-4*(m - 1))*R1_2 + (5*p*pow(m - 1, 2))*R1_3 + (-pow(p, 2)*pow(m - 1, 3))*R1_4;  // mu_phi_mu_phi_phi
            o[11] = (3*pow(m, 3))*R0_3 + (pow(m, 4)*p)*R0_4 + (3*pow(m - 1, 3))*R1_3 + (-p*pow(m - 1, 4))*R1_4;  // mu_phi_phi_phi_phi
            break;
        }
        case 4: {
            const double Q2_4 = R::psigamma(1.0 + X2, 4);
            const double R0_2 = d7::psi2_rest(X0);
            const double R0_3 = d7::psi3_rest(X0);
            const double R0_4 = d34_r4(X0);
            const double R1_2 = d7::psi2_rest(X1);
            const double R1_3 = d7::psi3_rest(X1);
            const double R1_4 = d34_r4(X1);
            o[0] = (-6*p*(2*m - 1)*(2*pow(m, 2) - 2*m + 1)/(pow(m, 4)*pow(m - 1, 4))) + (pow(p, 5))*R0_4 + (-pow(p, 5))*R1_4;  // mu_mu_mu_mu_mu
            o[1] = (2*(3*pow(m, 2) - 3*m + 1)/(pow(m, 3)*pow(m - 1, 3))) + (4*pow(p, 3))*R0_3 + (m*pow(p, 4))*R0_4 + (4*pow(p, 3))*R1_3 + (-pow(p, 4)*(m - 1))*R1_4;  // mu_mu_mu_mu_phi
            o[2] = (6*p)*R0_2 + (6*m*pow(p, 2))*R0_3 + (pow(m, 2)*pow(p, 3))*R0_4 + (-6*p)*R1_2 + (6*pow(p, 2)*(m - 1))*R1_3 + (-pow(p, 3)*pow(m - 1, 2))*R1_4;  // mu_mu_mu_phi_phi
            o[3] = (6*m)*R0_2 + (6*pow(m, 2)*p)*R0_3 + (pow(m, 3)*pow(p, 2))*R0_4 + (-6*(m - 1))*R1_2 + (6*p*pow(m - 1, 2))*R1_3 + (-pow(p, 2)*pow(m - 1, 3))*R1_4;  // mu_mu_phi_phi_phi
            o[4] = (6*p)*R0_2 + (6*m*pow(p, 2))*R0_3 + (pow(m, 2)*pow(p, 3))*R0_4 + (-6*p)*R1_2 + (6*pow(p, 2)*(m - 1))*R1_3 + (-pow(p, 3)*pow(m - 1, 2))*R1_4;  // phi_phi_mu_mu_mu
            o[5] = (6*m)*R0_2 + (6*pow(m, 2)*p)*R0_3 + (pow(m, 3)*pow(p, 2))*R0_4 + (-6*(m - 1))*R1_2 + (6*p*pow(m - 1, 2))*R1_3 + (-pow(p, 2)*pow(m - 1, 3))*R1_4;  // phi_phi_mu_mu_phi
            o[6] = (4*pow(m, 3))*R0_3 + (pow(m, 4)*p)*R0_4 + (4*pow(m - 1, 3))*R1_3 + (-p*pow(m - 1, 4))*R1_4;  // phi_phi_mu_phi_phi
            o[7] = (6*(p - 4)/pow(p, 5)) + (pow(m, 5))*R0_4 + (-pow(m - 1, 5))*R1_4 + (1)*Q2_4;  // phi_phi_phi_phi_phi
            o[8] = (3*pow(p, 3))*R0_3 + (m*pow(p, 4))*R0_4 + (3*pow(p, 3))*R1_3 + (-pow(p, 4)*(m - 1))*R1_4;  // mu_phi_mu_mu_mu
            o[9] = (4*p)*R0_2 + (5*m*pow(p, 2))*R0_3 + (pow(m, 2)*pow(p, 3))*R0_4 + (-4*p)*R1_2 + (5*pow(p, 2)*(m - 1))*R1_3 + (-pow(p, 3)*pow(m - 1, 2))*R1_4;  // mu_phi_mu_mu_phi
            o[10] = (4*m)*R0_2 + (5*pow(m, 2)*p)*R0_3 + (pow(m, 3)*pow(p, 2))*R0_4 + (-4*(m - 1))*R1_2 + (5*p*pow(m - 1, 2))*R1_3 + (-pow(p, 2)*pow(m - 1, 3))*R1_4;  // mu_phi_mu_phi_phi
            o[11] = (3*pow(m, 3))*R0_3 + (pow(m, 4)*p)*R0_4 + (3*pow(m - 1, 3))*R1_3 + (-p*pow(m - 1, 4))*R1_4;  // mu_phi_phi_phi_phi
            break;
        }
        case 2: {
            const double Q1_2 = R::psigamma(1.0 + X1, 2);
            const double Q1_3 = R::psigamma(1.0 + X1, 3);
            const double Q1_4 = R::psigamma(1.0 + X1, 4);
            const double R0_2 = d7::psi2_rest(X0);
            const double R0_3 = d7::psi3_rest(X0);
            const double R0_4 = d34_r4(X0);
            const double R2_4 = d34_r4(X2);
            o[0] = (6*(pow(m, 5)*p - 5*pow(m, 4)*p + 4*pow(m, 4) + 10*pow(m, 3)*p - 10*pow(m, 2)*p + 5*m*p - p)/(pow(m, 4)*pow(m - 1, 5))) + (pow(p, 5))*R0_4 + (pow(p, 5))*Q1_4;  // mu_mu_mu_mu_mu
            o[1] = (-2/pow(m, 3)) + (4*pow(p, 3))*R0_3 + (m*pow(p, 4))*R0_4 + (-4*pow(p, 3))*Q1_3 + (pow(p, 4)*(m - 1))*Q1_4;  // mu_mu_mu_mu_phi
            o[2] = (6*p)*R0_2 + (6*m*pow(p, 2))*R0_3 + (pow(m, 2)*pow(p, 3))*R0_4 + (6*p)*Q1_2 + (-6*pow(p, 2)*(m - 1))*Q1_3 + (pow(p, 3)*pow(m - 1, 2))*Q1_4;  // mu_mu_mu_phi_phi
            o[3] = (6*m)*R0_2 + (6*pow(m, 2)*p)*R0_3 + (pow(m, 3)*pow(p, 2))*R0_4 + (6*(m - 1))*Q1_2 + (-6*p*pow(m - 1, 2))*Q1_3 + (pow(p, 2)*pow(m - 1, 3))*Q1_4;  // mu_mu_phi_phi_phi
            o[4] = (6*p)*R0_2 + (6*m*pow(p, 2))*R0_3 + (pow(m, 2)*pow(p, 3))*R0_4 + (6*p)*Q1_2 + (-6*pow(p, 2)*(m - 1))*Q1_3 + (pow(p, 3)*pow(m - 1, 2))*Q1_4;  // phi_phi_mu_mu_mu
            o[5] = (6*m)*R0_2 + (6*pow(m, 2)*p)*R0_3 + (pow(m, 3)*pow(p, 2))*R0_4 + (6*(m - 1))*Q1_2 + (-6*p*pow(m - 1, 2))*Q1_3 + (pow(p, 2)*pow(m - 1, 3))*Q1_4;  // phi_phi_mu_mu_phi
            o[6] = (-2/pow(p, 3)) + (4*pow(m, 3))*R0_3 + (pow(m, 4)*p)*R0_4 + (-4*pow(m - 1, 3))*Q1_3 + (p*pow(m - 1, 4))*Q1_4;  // phi_phi_mu_phi_phi
            o[7] = (6*(m*p - p + 4)/pow(p, 5)) + (pow(m, 5))*R0_4 + (pow(m - 1, 5))*Q1_4 + (-1)*R2_4;  // phi_phi_phi_phi_phi
            o[8] = (6/(p*pow(m - 1, 4))) + (3*pow(p, 3))*R0_3 + (m*pow(p, 4))*R0_4 + (-3*pow(p, 3))*Q1_3 + (pow(p, 4)*(m - 1))*Q1_4;  // mu_phi_mu_mu_mu
            o[9] = (2/(pow(p, 2)*pow(m - 1, 3))) + (4*p)*R0_2 + (5*m*pow(p, 2))*R0_3 + (pow(m, 2)*pow(p, 3))*R0_4 + (4*p)*Q1_2 + (-5*pow(p, 2)*(m - 1))*Q1_3 + (pow(p, 3)*pow(m - 1, 2))*Q1_4;  // mu_phi_mu_mu_phi
            o[10] = (2/(pow(p, 3)*pow(m - 1, 2))) + (4*m)*R0_2 + (5*pow(m, 2)*p)*R0_3 + (pow(m, 3)*pow(p, 2))*R0_4 + (4*(m - 1))*Q1_2 + (-5*p*pow(m - 1, 2))*Q1_3 + (pow(p, 2)*pow(m - 1, 3))*Q1_4;  // mu_phi_mu_phi_phi
            o[11] = (6/(pow(p, 4)*(m - 1))) + (3*pow(m, 3))*R0_3 + (pow(m, 4)*p)*R0_4 + (-3*pow(m - 1, 3))*Q1_3 + (p*pow(m - 1, 4))*Q1_4;  // mu_phi_phi_phi_phi
            break;
        }
        case 6: {
            const double Q1_2 = R::psigamma(1.0 + X1, 2);
            const double Q1_3 = R::psigamma(1.0 + X1, 3);
            const double Q1_4 = R::psigamma(1.0 + X1, 4);
            const double Q2_4 = R::psigamma(1.0 + X2, 4);
            const double R0_2 = d7::psi2_rest(X0);
            const double R0_3 = d7::psi3_rest(X0);
            const double R0_4 = d34_r4(X0);
            o[0] = (6*(pow(m, 5)*p - 5*pow(m, 4)*p + 4*pow(m, 4) + 10*pow(m, 3)*p - 10*pow(m, 2)*p + 5*m*p - p)/(pow(m, 4)*pow(m - 1, 5))) + (pow(p, 5))*R0_4 + (pow(p, 5))*Q1_4;  // mu_mu_mu_mu_mu
            o[1] = (-2/pow(m, 3)) + (4*pow(p, 3))*R0_3 + (m*pow(p, 4))*R0_4 + (-4*pow(p, 3))*Q1_3 + (pow(p, 4)*(m - 1))*Q1_4;  // mu_mu_mu_mu_phi
            o[2] = (6*p)*R0_2 + (6*m*pow(p, 2))*R0_3 + (pow(m, 2)*pow(p, 3))*R0_4 + (6*p)*Q1_2 + (-6*pow(p, 2)*(m - 1))*Q1_3 + (pow(p, 3)*pow(m - 1, 2))*Q1_4;  // mu_mu_mu_phi_phi
            o[3] = (6*m)*R0_2 + (6*pow(m, 2)*p)*R0_3 + (pow(m, 3)*pow(p, 2))*R0_4 + (6*(m - 1))*Q1_2 + (-6*p*pow(m - 1, 2))*Q1_3 + (pow(p, 2)*pow(m - 1, 3))*Q1_4;  // mu_mu_phi_phi_phi
            o[4] = (6*p)*R0_2 + (6*m*pow(p, 2))*R0_3 + (pow(m, 2)*pow(p, 3))*R0_4 + (6*p)*Q1_2 + (-6*pow(p, 2)*(m - 1))*Q1_3 + (pow(p, 3)*pow(m - 1, 2))*Q1_4;  // phi_phi_mu_mu_mu
            o[5] = (6*m)*R0_2 + (6*pow(m, 2)*p)*R0_3 + (pow(m, 3)*pow(p, 2))*R0_4 + (6*(m - 1))*Q1_2 + (-6*p*pow(m - 1, 2))*Q1_3 + (pow(p, 2)*pow(m - 1, 3))*Q1_4;  // phi_phi_mu_mu_phi
            o[6] = (-2/pow(p, 3)) + (4*pow(m, 3))*R0_3 + (pow(m, 4)*p)*R0_4 + (-4*pow(m - 1, 3))*Q1_3 + (p*pow(m - 1, 4))*Q1_4;  // phi_phi_mu_phi_phi
            o[7] = (6*m/pow(p, 4)) + (pow(m, 5))*R0_4 + (pow(m - 1, 5))*Q1_4 + (1)*Q2_4;  // phi_phi_phi_phi_phi
            o[8] = (6/(p*pow(m - 1, 4))) + (3*pow(p, 3))*R0_3 + (m*pow(p, 4))*R0_4 + (-3*pow(p, 3))*Q1_3 + (pow(p, 4)*(m - 1))*Q1_4;  // mu_phi_mu_mu_mu
            o[9] = (2/(pow(p, 2)*pow(m - 1, 3))) + (4*p)*R0_2 + (5*m*pow(p, 2))*R0_3 + (pow(m, 2)*pow(p, 3))*R0_4 + (4*p)*Q1_2 + (-5*pow(p, 2)*(m - 1))*Q1_3 + (pow(p, 3)*pow(m - 1, 2))*Q1_4;  // mu_phi_mu_mu_phi
            o[10] = (2/(pow(p, 3)*pow(m - 1, 2))) + (4*m)*R0_2 + (5*pow(m, 2)*p)*R0_3 + (pow(m, 3)*pow(p, 2))*R0_4 + (4*(m - 1))*Q1_2 + (-5*p*pow(m - 1, 2))*Q1_3 + (pow(p, 2)*pow(m - 1, 3))*Q1_4;  // mu_phi_mu_phi_phi
            o[11] = (6/(pow(p, 4)*(m - 1))) + (3*pow(m, 3))*R0_3 + (pow(m, 4)*p)*R0_4 + (-3*pow(m - 1, 3))*Q1_3 + (p*pow(m - 1, 4))*Q1_4;  // mu_phi_phi_phi_phi
            break;
        }
        case 1: {
            const double Q0_2 = R::psigamma(1.0 + X0, 2);
            const double Q0_3 = R::psigamma(1.0 + X0, 3);
            const double Q0_4 = R::psigamma(1.0 + X0, 4);
            const double R1_2 = d7::psi2_rest(X1);
            const double R1_3 = d7::psi3_rest(X1);
            const double R1_4 = d34_r4(X1);
            const double R2_4 = d34_r4(X2);
            o[0] = (-6*(pow(m, 5)*p - 4*pow(m, 4) + 16*pow(m, 3) - 24*pow(m, 2) + 16*m - 4)/(pow(m, 5)*pow(m - 1, 4))) + (-pow(p, 5))*Q0_4 + (-pow(p, 5))*R1_4;  // mu_mu_mu_mu_mu
            o[1] = (2/pow(m - 1, 3)) + (-4*pow(p, 3))*Q0_3 + (-m*pow(p, 4))*Q0_4 + (4*pow(p, 3))*R1_3 + (-pow(p, 4)*(m - 1))*R1_4;  // mu_mu_mu_mu_phi
            o[2] = (-6*p)*Q0_2 + (-6*m*pow(p, 2))*Q0_3 + (-pow(m, 2)*pow(p, 3))*Q0_4 + (-6*p)*R1_2 + (6*pow(p, 2)*(m - 1))*R1_3 + (-pow(p, 3)*pow(m - 1, 2))*R1_4;  // mu_mu_mu_phi_phi
            o[3] = (-6*m)*Q0_2 + (-6*pow(m, 2)*p)*Q0_3 + (-pow(m, 3)*pow(p, 2))*Q0_4 + (-6*(m - 1))*R1_2 + (6*p*pow(m - 1, 2))*R1_3 + (-pow(p, 2)*pow(m - 1, 3))*R1_4;  // mu_mu_phi_phi_phi
            o[4] = (-6*p)*Q0_2 + (-6*m*pow(p, 2))*Q0_3 + (-pow(m, 2)*pow(p, 3))*Q0_4 + (-6*p)*R1_2 + (6*pow(p, 2)*(m - 1))*R1_3 + (-pow(p, 3)*pow(m - 1, 2))*R1_4;  // phi_phi_mu_mu_mu
            o[5] = (-6*m)*Q0_2 + (-6*pow(m, 2)*p)*Q0_3 + (-pow(m, 3)*pow(p, 2))*Q0_4 + (-6*(m - 1))*R1_2 + (6*p*pow(m - 1, 2))*R1_3 + (-pow(p, 2)*pow(m - 1, 3))*R1_4;  // phi_phi_mu_mu_phi
            o[6] = (2/pow(p, 3)) + (-4*pow(m, 3))*Q0_3 + (-pow(m, 4)*p)*Q0_4 + (4*pow(m - 1, 3))*R1_3 + (-p*pow(m - 1, 4))*R1_4;  // phi_phi_mu_phi_phi
            o[7] = (-6*(m*p - 4)/pow(p, 5)) + (-pow(m, 5))*Q0_4 + (-pow(m - 1, 5))*R1_4 + (-1)*R2_4;  // phi_phi_phi_phi_phi
            o[8] = (6/(pow(m, 4)*p)) + (-3*pow(p, 3))*Q0_3 + (-m*pow(p, 4))*Q0_4 + (3*pow(p, 3))*R1_3 + (-pow(p, 4)*(m - 1))*R1_4;  // mu_phi_mu_mu_mu
            o[9] = (2/(pow(m, 3)*pow(p, 2))) + (-4*p)*Q0_2 + (-5*m*pow(p, 2))*Q0_3 + (-pow(m, 2)*pow(p, 3))*Q0_4 + (-4*p)*R1_2 + (5*pow(p, 2)*(m - 1))*R1_3 + (-pow(p, 3)*pow(m - 1, 2))*R1_4;  // mu_phi_mu_mu_phi
            o[10] = (2/(pow(m, 2)*pow(p, 3))) + (-4*m)*Q0_2 + (-5*pow(m, 2)*p)*Q0_3 + (-pow(m, 3)*pow(p, 2))*Q0_4 + (-4*(m - 1))*R1_2 + (5*p*pow(m - 1, 2))*R1_3 + (-pow(p, 2)*pow(m - 1, 3))*R1_4;  // mu_phi_mu_phi_phi
            o[11] = (6/(m*pow(p, 4))) + (-3*pow(m, 3))*Q0_3 + (-pow(m, 4)*p)*Q0_4 + (3*pow(m - 1, 3))*R1_3 + (-p*pow(m - 1, 4))*R1_4;  // mu_phi_phi_phi_phi
            break;
        }
        case 5: {
            const double Q0_2 = R::psigamma(1.0 + X0, 2);
            const double Q0_3 = R::psigamma(1.0 + X0, 3);
            const double Q0_4 = R::psigamma(1.0 + X0, 4);
            const double Q2_4 = R::psigamma(1.0 + X2, 4);
            const double R1_2 = d7::psi2_rest(X1);
            const double R1_3 = d7::psi3_rest(X1);
            const double R1_4 = d34_r4(X1);
            o[0] = (-6*(pow(m, 5)*p - 4*pow(m, 4) + 16*pow(m, 3) - 24*pow(m, 2) + 16*m - 4)/(pow(m, 5)*pow(m - 1, 4))) + (-pow(p, 5))*Q0_4 + (-pow(p, 5))*R1_4;  // mu_mu_mu_mu_mu
            o[1] = (2/pow(m - 1, 3)) + (-4*pow(p, 3))*Q0_3 + (-m*pow(p, 4))*Q0_4 + (4*pow(p, 3))*R1_3 + (-pow(p, 4)*(m - 1))*R1_4;  // mu_mu_mu_mu_phi
            o[2] = (-6*p)*Q0_2 + (-6*m*pow(p, 2))*Q0_3 + (-pow(m, 2)*pow(p, 3))*Q0_4 + (-6*p)*R1_2 + (6*pow(p, 2)*(m - 1))*R1_3 + (-pow(p, 3)*pow(m - 1, 2))*R1_4;  // mu_mu_mu_phi_phi
            o[3] = (-6*m)*Q0_2 + (-6*pow(m, 2)*p)*Q0_3 + (-pow(m, 3)*pow(p, 2))*Q0_4 + (-6*(m - 1))*R1_2 + (6*p*pow(m - 1, 2))*R1_3 + (-pow(p, 2)*pow(m - 1, 3))*R1_4;  // mu_mu_phi_phi_phi
            o[4] = (-6*p)*Q0_2 + (-6*m*pow(p, 2))*Q0_3 + (-pow(m, 2)*pow(p, 3))*Q0_4 + (-6*p)*R1_2 + (6*pow(p, 2)*(m - 1))*R1_3 + (-pow(p, 3)*pow(m - 1, 2))*R1_4;  // phi_phi_mu_mu_mu
            o[5] = (-6*m)*Q0_2 + (-6*pow(m, 2)*p)*Q0_3 + (-pow(m, 3)*pow(p, 2))*Q0_4 + (-6*(m - 1))*R1_2 + (6*p*pow(m - 1, 2))*R1_3 + (-pow(p, 2)*pow(m - 1, 3))*R1_4;  // phi_phi_mu_mu_phi
            o[6] = (2/pow(p, 3)) + (-4*pow(m, 3))*Q0_3 + (-pow(m, 4)*p)*Q0_4 + (4*pow(m - 1, 3))*R1_3 + (-p*pow(m - 1, 4))*R1_4;  // phi_phi_mu_phi_phi
            o[7] = (-6*(m - 1)/pow(p, 4)) + (-pow(m, 5))*Q0_4 + (-pow(m - 1, 5))*R1_4 + (1)*Q2_4;  // phi_phi_phi_phi_phi
            o[8] = (6/(pow(m, 4)*p)) + (-3*pow(p, 3))*Q0_3 + (-m*pow(p, 4))*Q0_4 + (3*pow(p, 3))*R1_3 + (-pow(p, 4)*(m - 1))*R1_4;  // mu_phi_mu_mu_mu
            o[9] = (2/(pow(m, 3)*pow(p, 2))) + (-4*p)*Q0_2 + (-5*m*pow(p, 2))*Q0_3 + (-pow(m, 2)*pow(p, 3))*Q0_4 + (-4*p)*R1_2 + (5*pow(p, 2)*(m - 1))*R1_3 + (-pow(p, 3)*pow(m - 1, 2))*R1_4;  // mu_phi_mu_mu_phi
            o[10] = (2/(pow(m, 2)*pow(p, 3))) + (-4*m)*Q0_2 + (-5*pow(m, 2)*p)*Q0_3 + (-pow(m, 3)*pow(p, 2))*Q0_4 + (-4*(m - 1))*R1_2 + (5*p*pow(m - 1, 2))*R1_3 + (-pow(p, 2)*pow(m - 1, 3))*R1_4;  // mu_phi_mu_phi_phi
            o[11] = (6/(m*pow(p, 4))) + (-3*pow(m, 3))*Q0_3 + (-pow(m, 4)*p)*Q0_4 + (3*pow(m - 1, 3))*R1_3 + (-p*pow(m - 1, 4))*R1_4;  // mu_phi_phi_phi_phi
            break;
        }
        case 3: {
            const double Q0_2 = R::psigamma(1.0 + X0, 2);
            const double Q0_3 = R::psigamma(1.0 + X0, 3);
            const double Q0_4 = R::psigamma(1.0 + X0, 4);
            const double Q1_2 = R::psigamma(1.0 + X1, 2);
            const double Q1_3 = R::psigamma(1.0 + X1, 3);
            const double Q1_4 = R::psigamma(1.0 + X1, 4);
            const double R2_4 = d34_r4(X2);
            o[0] = (24*(2*m - 1)*(pow(m, 4) - 2*pow(m, 3) + 4*pow(m, 2) - 3*m + 1)/(pow(m, 5)*pow(m - 1, 5))) + (-pow(p, 5))*Q0_4 + (pow(p, 5))*Q1_4;  // mu_mu_mu_mu_mu
            o[1] = (-4*pow(p, 3))*Q0_3 + (-m*pow(p, 4))*Q0_4 + (-4*pow(p, 3))*Q1_3 + (pow(p, 4)*(m - 1))*Q1_4;  // mu_mu_mu_mu_phi
            o[2] = (-6*p)*Q0_2 + (-6*m*pow(p, 2))*Q0_3 + (-pow(m, 2)*pow(p, 3))*Q0_4 + (6*p)*Q1_2 + (-6*pow(p, 2)*(m - 1))*Q1_3 + (pow(p, 3)*pow(m - 1, 2))*Q1_4;  // mu_mu_mu_phi_phi
            o[3] = (-6*m)*Q0_2 + (-6*pow(m, 2)*p)*Q0_3 + (-pow(m, 3)*pow(p, 2))*Q0_4 + (6*(m - 1))*Q1_2 + (-6*p*pow(m - 1, 2))*Q1_3 + (pow(p, 2)*pow(m - 1, 3))*Q1_4;  // mu_mu_phi_phi_phi
            o[4] = (-6*p)*Q0_2 + (-6*m*pow(p, 2))*Q0_3 + (-pow(m, 2)*pow(p, 3))*Q0_4 + (6*p)*Q1_2 + (-6*pow(p, 2)*(m - 1))*Q1_3 + (pow(p, 3)*pow(m - 1, 2))*Q1_4;  // phi_phi_mu_mu_mu
            o[5] = (-6*m)*Q0_2 + (-6*pow(m, 2)*p)*Q0_3 + (-pow(m, 3)*pow(p, 2))*Q0_4 + (6*(m - 1))*Q1_2 + (-6*p*pow(m - 1, 2))*Q1_3 + (pow(p, 2)*pow(m - 1, 3))*Q1_4;  // phi_phi_mu_mu_phi
            o[6] = (-4*pow(m, 3))*Q0_3 + (-pow(m, 4)*p)*Q0_4 + (-4*pow(m - 1, 3))*Q1_3 + (p*pow(m - 1, 4))*Q1_4;  // phi_phi_mu_phi_phi
            o[7] = (-6*(p - 8)/pow(p, 5)) + (-pow(m, 5))*Q0_4 + (pow(m - 1, 5))*Q1_4 + (-1)*R2_4;  // phi_phi_phi_phi_phi
            o[8] = (6*(2*pow(m, 4) - 4*pow(m, 3) + 6*pow(m, 2) - 4*m + 1)/(pow(m, 4)*p*pow(m - 1, 4))) + (-3*pow(p, 3))*Q0_3 + (-m*pow(p, 4))*Q0_4 + (-3*pow(p, 3))*Q1_3 + (pow(p, 4)*(m - 1))*Q1_4;  // mu_phi_mu_mu_mu
            o[9] = (2*(2*m - 1)*(pow(m, 2) - m + 1)/(pow(m, 3)*pow(p, 2)*pow(m - 1, 3))) + (-4*p)*Q0_2 + (-5*m*pow(p, 2))*Q0_3 + (-pow(m, 2)*pow(p, 3))*Q0_4 + (4*p)*Q1_2 + (-5*pow(p, 2)*(m - 1))*Q1_3 + (pow(p, 3)*pow(m - 1, 2))*Q1_4;  // mu_phi_mu_mu_phi
            o[10] = (2*(2*pow(m, 2) - 2*m + 1)/(pow(m, 2)*pow(p, 3)*pow(m - 1, 2))) + (-4*m)*Q0_2 + (-5*pow(m, 2)*p)*Q0_3 + (-pow(m, 3)*pow(p, 2))*Q0_4 + (4*(m - 1))*Q1_2 + (-5*p*pow(m - 1, 2))*Q1_3 + (pow(p, 2)*pow(m - 1, 3))*Q1_4;  // mu_phi_mu_phi_phi
            o[11] = (6*(2*m - 1)/(m*pow(p, 4)*(m - 1))) + (-3*pow(m, 3))*Q0_3 + (-pow(m, 4)*p)*Q0_4 + (-3*pow(m - 1, 3))*Q1_3 + (p*pow(m - 1, 4))*Q1_4;  // mu_phi_phi_phi_phi
            break;
        }
        case 7: {
            const double Q0_2 = R::psigamma(1.0 + X0, 2);
            const double Q0_3 = R::psigamma(1.0 + X0, 3);
            const double Q0_4 = R::psigamma(1.0 + X0, 4);
            const double Q1_2 = R::psigamma(1.0 + X1, 2);
            const double Q1_3 = R::psigamma(1.0 + X1, 3);
            const double Q1_4 = R::psigamma(1.0 + X1, 4);
            const double Q2_4 = R::psigamma(1.0 + X2, 4);
            o[0] = (24*(2*m - 1)*(pow(m, 4) - 2*pow(m, 3) + 4*pow(m, 2) - 3*m + 1)/(pow(m, 5)*pow(m - 1, 5))) + (-pow(p, 5))*Q0_4 + (pow(p, 5))*Q1_4;  // mu_mu_mu_mu_mu
            o[1] = (-4*pow(p, 3))*Q0_3 + (-m*pow(p, 4))*Q0_4 + (-4*pow(p, 3))*Q1_3 + (pow(p, 4)*(m - 1))*Q1_4;  // mu_mu_mu_mu_phi
            o[2] = (-6*p)*Q0_2 + (-6*m*pow(p, 2))*Q0_3 + (-pow(m, 2)*pow(p, 3))*Q0_4 + (6*p)*Q1_2 + (-6*pow(p, 2)*(m - 1))*Q1_3 + (pow(p, 3)*pow(m - 1, 2))*Q1_4;  // mu_mu_mu_phi_phi
            o[3] = (-6*m)*Q0_2 + (-6*pow(m, 2)*p)*Q0_3 + (-pow(m, 3)*pow(p, 2))*Q0_4 + (6*(m - 1))*Q1_2 + (-6*p*pow(m - 1, 2))*Q1_3 + (pow(p, 2)*pow(m - 1, 3))*Q1_4;  // mu_mu_phi_phi_phi
            o[4] = (-6*p)*Q0_2 + (-6*m*pow(p, 2))*Q0_3 + (-pow(m, 2)*pow(p, 3))*Q0_4 + (6*p)*Q1_2 + (-6*pow(p, 2)*(m - 1))*Q1_3 + (pow(p, 3)*pow(m - 1, 2))*Q1_4;  // phi_phi_mu_mu_mu
            o[5] = (-6*m)*Q0_2 + (-6*pow(m, 2)*p)*Q0_3 + (-pow(m, 3)*pow(p, 2))*Q0_4 + (6*(m - 1))*Q1_2 + (-6*p*pow(m - 1, 2))*Q1_3 + (pow(p, 2)*pow(m - 1, 3))*Q1_4;  // phi_phi_mu_mu_phi
            o[6] = (-4*pow(m, 3))*Q0_3 + (-pow(m, 4)*p)*Q0_4 + (-4*pow(m - 1, 3))*Q1_3 + (p*pow(m - 1, 4))*Q1_4;  // phi_phi_mu_phi_phi
            o[7] = (24/pow(p, 5)) + (-pow(m, 5))*Q0_4 + (pow(m - 1, 5))*Q1_4 + (1)*Q2_4;  // phi_phi_phi_phi_phi
            o[8] = (6*(2*pow(m, 4) - 4*pow(m, 3) + 6*pow(m, 2) - 4*m + 1)/(pow(m, 4)*p*pow(m - 1, 4))) + (-3*pow(p, 3))*Q0_3 + (-m*pow(p, 4))*Q0_4 + (-3*pow(p, 3))*Q1_3 + (pow(p, 4)*(m - 1))*Q1_4;  // mu_phi_mu_mu_mu
            o[9] = (2*(2*m - 1)*(pow(m, 2) - m + 1)/(pow(m, 3)*pow(p, 2)*pow(m - 1, 3))) + (-4*p)*Q0_2 + (-5*m*pow(p, 2))*Q0_3 + (-pow(m, 2)*pow(p, 3))*Q0_4 + (4*p)*Q1_2 + (-5*pow(p, 2)*(m - 1))*Q1_3 + (pow(p, 3)*pow(m - 1, 2))*Q1_4;  // mu_phi_mu_mu_phi
            o[10] = (2*(2*pow(m, 2) - 2*m + 1)/(pow(m, 2)*pow(p, 3)*pow(m - 1, 2))) + (-4*m)*Q0_2 + (-5*pow(m, 2)*p)*Q0_3 + (-pow(m, 3)*pow(p, 2))*Q0_4 + (4*(m - 1))*Q1_2 + (-5*p*pow(m - 1, 2))*Q1_3 + (pow(p, 2)*pow(m - 1, 3))*Q1_4;  // mu_phi_mu_phi_phi
            o[11] = (6*(2*m - 1)/(m*pow(p, 4)*(m - 1))) + (-3*pow(m, 3))*Q0_3 + (-pow(m, 4)*p)*Q0_4 + (-3*pow(m - 1, 3))*Q1_3 + (p*pow(m - 1, 4))*Q1_4;  // mu_phi_phi_phi_phi
            break;
        }
        }
    });
}

// [[Rcpp::export]]
List beta1_dexpected4_cpp(NumericVector y, NumericVector mu, NumericVector phi, int threads = 1) {
    d34e::Par A0(mu), A1(phi);
    return d34e::run(y.size(), threads, d7::kMinCostly, 15,
                     [&](std::size_t i, d34e::Row o) {
        const double m = A0[i], p = A1[i];
        const double X0 = m*p;
        const double X1 = p*(1 - m);
        const double X2 = p;
        switch ((X0 < 1.0 ? 1 : 0) | (X1 < 1.0 ? 2 : 0) | (X2 < 1.0 ? 4 : 0)) {
        case 0: {
            const double R0_2 = d7::psi2_rest(X0);
            const double R0_3 = d7::psi3_rest(X0);
            const double R0_4 = d34_r4(X0);
            const double R0_5 = d34_r5(X0);
            const double R1_2 = d7::psi2_rest(X1);
            const double R1_3 = d7::psi3_rest(X1);
            const double R1_4 = d34_r4(X1);
            const double R1_5 = d34_r5(X1);
            const double R2_5 = d34_r5(X2);
            o[0] = (24*p*(5*pow(m, 4) - 10*pow(m, 3) + 10*pow(m, 2) - 5*m + 1)/(pow(m, 5)*pow(m - 1, 5))) + (pow(p, 6))*R0_5 + (pow(p, 6))*R1_5;  // mu_mu_mu_mu_mu_mu
            o[1] = (-6*(2*m - 1)*(2*pow(m, 2) - 2*m + 1)/(pow(m, 4)*pow(m - 1, 4))) + (5*pow(p, 4))*R0_4 + (m*pow(p, 5))*R0_5 + (-5*pow(p, 4))*R1_4 + (pow(p, 5)*(m - 1))*R1_5;  // mu_mu_mu_mu_mu_phi
            o[2] = (12*pow(p, 2))*R0_3 + (8*m*pow(p, 3))*R0_4 + (pow(m, 2)*pow(p, 4))*R0_5 + (12*pow(p, 2))*R1_3 + (-8*pow(p, 3)*(m - 1))*R1_4 + (pow(p, 4)*pow(m - 1, 2))*R1_5;  // mu_mu_mu_mu_phi_phi
            o[3] = (6)*R0_2 + (18*m*p)*R0_3 + (9*pow(m, 2)*pow(p, 2))*R0_4 + (pow(m, 3)*pow(p, 3))*R0_5 + (-6)*R1_2 + (18*p*(m - 1))*R1_3 + (-9*pow(p, 2)*pow(m - 1, 2))*R1_4 + (pow(p, 3)*pow(m - 1, 3))*R1_5;  // mu_mu_mu_phi_phi_phi
            o[4] = (12*pow(m, 2))*R0_3 + (8*pow(m, 3)*p)*R0_4 + (pow(m, 4)*pow(p, 2))*R0_5 + (12*pow(m - 1, 2))*R1_3 + (-8*p*pow(m - 1, 3))*R1_4 + (pow(p, 2)*pow(m - 1, 4))*R1_5;  // mu_mu_phi_phi_phi_phi
            o[5] = (12*pow(p, 2))*R0_3 + (8*m*pow(p, 3))*R0_4 + (pow(m, 2)*pow(p, 4))*R0_5 + (12*pow(p, 2))*R1_3 + (-8*pow(p, 3)*(m - 1))*R1_4 + (pow(p, 4)*pow(m - 1, 2))*R1_5;  // phi_phi_mu_mu_mu_mu
            o[6] = (6)*R0_2 + (18*m*p)*R0_3 + (9*pow(m, 2)*pow(p, 2))*R0_4 + (pow(m, 3)*pow(p, 3))*R0_5 + (-6)*R1_2 + (18*p*(m - 1))*R1_3 + (-9*pow(p, 2)*pow(m - 1, 2))*R1_4 + (pow(p, 3)*pow(m - 1, 3))*R1_5;  // phi_phi_mu_mu_mu_phi
            o[7] = (12*pow(m, 2))*R0_3 + (8*pow(m, 3)*p)*R0_4 + (pow(m, 4)*pow(p, 2))*R0_5 + (12*pow(m - 1, 2))*R1_3 + (-8*p*pow(m - 1, 3))*R1_4 + (pow(p, 2)*pow(m - 1, 4))*R1_5;  // phi_phi_mu_mu_phi_phi
            o[8] = (5*pow(m, 4))*R0_4 + (pow(m, 5)*p)*R0_5 + (-5*pow(m - 1, 4))*R1_4 + (p*pow(m - 1, 5))*R1_5;  // phi_phi_mu_phi_phi_phi
            o[9] = (pow(m, 6))*R0_5 + (pow(m - 1, 6))*R1_5 + (-1)*R2_5;  // phi_phi_phi_phi_phi_phi
            o[10] = (4*pow(p, 4))*R0_4 + (m*pow(p, 5))*R0_5 + (-4*pow(p, 4))*R1_4 + (pow(p, 5)*(m - 1))*R1_5;  // mu_phi_mu_mu_mu_mu
            o[11] = (9*pow(p, 2))*R0_3 + (7*m*pow(p, 3))*R0_4 + (pow(m, 2)*pow(p, 4))*R0_5 + (9*pow(p, 2))*R1_3 + (-7*pow(p, 3)*(m - 1))*R1_4 + (pow(p, 4)*pow(m - 1, 2))*R1_5;  // mu_phi_mu_mu_mu_phi
            o[12] = (4)*R0_2 + (14*m*p)*R0_3 + (8*pow(m, 2)*pow(p, 2))*R0_4 + (pow(m, 3)*pow(p, 3))*R0_5 + (-4)*R1_2 + (14*p*(m - 1))*R1_3 + (-8*pow(p, 2)*pow(m - 1, 2))*R1_4 + (pow(p, 3)*pow(m - 1, 3))*R1_5;  // mu_phi_mu_mu_phi_phi
            o[13] = (9*pow(m, 2))*R0_3 + (7*pow(m, 3)*p)*R0_4 + (pow(m, 4)*pow(p, 2))*R0_5 + (9*pow(m - 1, 2))*R1_3 + (-7*p*pow(m - 1, 3))*R1_4 + (pow(p, 2)*pow(m - 1, 4))*R1_5;  // mu_phi_mu_phi_phi_phi
            o[14] = (4*pow(m, 4))*R0_4 + (pow(m, 5)*p)*R0_5 + (-4*pow(m - 1, 4))*R1_4 + (p*pow(m - 1, 5))*R1_5;  // mu_phi_phi_phi_phi_phi
            break;
        }
        case 4: {
            const double Q2_5 = R::psigamma(1.0 + X2, 5);
            const double R0_2 = d7::psi2_rest(X0);
            const double R0_3 = d7::psi3_rest(X0);
            const double R0_4 = d34_r4(X0);
            const double R0_5 = d34_r5(X0);
            const double R1_2 = d7::psi2_rest(X1);
            const double R1_3 = d7::psi3_rest(X1);
            const double R1_4 = d34_r4(X1);
            const double R1_5 = d34_r5(X1);
            o[0] = (24*p*(5*pow(m, 4) - 10*pow(m, 3) + 10*pow(m, 2) - 5*m + 1)/(pow(m, 5)*pow(m - 1, 5))) + (pow(p, 6))*R0_5 + (pow(p, 6))*R1_5;  // mu_mu_mu_mu_mu_mu
            o[1] = (-6*(2*m - 1)*(2*pow(m, 2) - 2*m + 1)/(pow(m, 4)*pow(m - 1, 4))) + (5*pow(p, 4))*R0_4 + (m*pow(p, 5))*R0_5 + (-5*pow(p, 4))*R1_4 + (pow(p, 5)*(m - 1))*R1_5;  // mu_mu_mu_mu_mu_phi
            o[2] = (12*pow(p, 2))*R0_3 + (8*m*pow(p, 3))*R0_4 + (pow(m, 2)*pow(p, 4))*R0_5 + (12*pow(p, 2))*R1_3 + (-8*pow(p, 3)*(m - 1))*R1_4 + (pow(p, 4)*pow(m - 1, 2))*R1_5;  // mu_mu_mu_mu_phi_phi
            o[3] = (6)*R0_2 + (18*m*p)*R0_3 + (9*pow(m, 2)*pow(p, 2))*R0_4 + (pow(m, 3)*pow(p, 3))*R0_5 + (-6)*R1_2 + (18*p*(m - 1))*R1_3 + (-9*pow(p, 2)*pow(m - 1, 2))*R1_4 + (pow(p, 3)*pow(m - 1, 3))*R1_5;  // mu_mu_mu_phi_phi_phi
            o[4] = (12*pow(m, 2))*R0_3 + (8*pow(m, 3)*p)*R0_4 + (pow(m, 4)*pow(p, 2))*R0_5 + (12*pow(m - 1, 2))*R1_3 + (-8*p*pow(m - 1, 3))*R1_4 + (pow(p, 2)*pow(m - 1, 4))*R1_5;  // mu_mu_phi_phi_phi_phi
            o[5] = (12*pow(p, 2))*R0_3 + (8*m*pow(p, 3))*R0_4 + (pow(m, 2)*pow(p, 4))*R0_5 + (12*pow(p, 2))*R1_3 + (-8*pow(p, 3)*(m - 1))*R1_4 + (pow(p, 4)*pow(m - 1, 2))*R1_5;  // phi_phi_mu_mu_mu_mu
            o[6] = (6)*R0_2 + (18*m*p)*R0_3 + (9*pow(m, 2)*pow(p, 2))*R0_4 + (pow(m, 3)*pow(p, 3))*R0_5 + (-6)*R1_2 + (18*p*(m - 1))*R1_3 + (-9*pow(p, 2)*pow(m - 1, 2))*R1_4 + (pow(p, 3)*pow(m - 1, 3))*R1_5;  // phi_phi_mu_mu_mu_phi
            o[7] = (12*pow(m, 2))*R0_3 + (8*pow(m, 3)*p)*R0_4 + (pow(m, 4)*pow(p, 2))*R0_5 + (12*pow(m - 1, 2))*R1_3 + (-8*p*pow(m - 1, 3))*R1_4 + (pow(p, 2)*pow(m - 1, 4))*R1_5;  // phi_phi_mu_mu_phi_phi
            o[8] = (5*pow(m, 4))*R0_4 + (pow(m, 5)*p)*R0_5 + (-5*pow(m - 1, 4))*R1_4 + (p*pow(m - 1, 5))*R1_5;  // phi_phi_mu_phi_phi_phi
            o[9] = (-24*(p - 5)/pow(p, 6)) + (pow(m, 6))*R0_5 + (pow(m - 1, 6))*R1_5 + (1)*Q2_5;  // phi_phi_phi_phi_phi_phi
            o[10] = (4*pow(p, 4))*R0_4 + (m*pow(p, 5))*R0_5 + (-4*pow(p, 4))*R1_4 + (pow(p, 5)*(m - 1))*R1_5;  // mu_phi_mu_mu_mu_mu
            o[11] = (9*pow(p, 2))*R0_3 + (7*m*pow(p, 3))*R0_4 + (pow(m, 2)*pow(p, 4))*R0_5 + (9*pow(p, 2))*R1_3 + (-7*pow(p, 3)*(m - 1))*R1_4 + (pow(p, 4)*pow(m - 1, 2))*R1_5;  // mu_phi_mu_mu_mu_phi
            o[12] = (4)*R0_2 + (14*m*p)*R0_3 + (8*pow(m, 2)*pow(p, 2))*R0_4 + (pow(m, 3)*pow(p, 3))*R0_5 + (-4)*R1_2 + (14*p*(m - 1))*R1_3 + (-8*pow(p, 2)*pow(m - 1, 2))*R1_4 + (pow(p, 3)*pow(m - 1, 3))*R1_5;  // mu_phi_mu_mu_phi_phi
            o[13] = (9*pow(m, 2))*R0_3 + (7*pow(m, 3)*p)*R0_4 + (pow(m, 4)*pow(p, 2))*R0_5 + (9*pow(m - 1, 2))*R1_3 + (-7*p*pow(m - 1, 3))*R1_4 + (pow(p, 2)*pow(m - 1, 4))*R1_5;  // mu_phi_mu_phi_phi_phi
            o[14] = (4*pow(m, 4))*R0_4 + (pow(m, 5)*p)*R0_5 + (-4*pow(m - 1, 4))*R1_4 + (p*pow(m - 1, 5))*R1_5;  // mu_phi_phi_phi_phi_phi
            break;
        }
        case 2: {
            const double Q1_2 = R::psigamma(1.0 + X1, 2);
            const double Q1_3 = R::psigamma(1.0 + X1, 3);
            const double Q1_4 = R::psigamma(1.0 + X1, 4);
            const double Q1_5 = R::psigamma(1.0 + X1, 5);
            const double R0_2 = d7::psi2_rest(X0);
            const double R0_3 = d7::psi3_rest(X0);
            const double R0_4 = d34_r4(X0);
            const double R0_5 = d34_r5(X0);
            const double R2_5 = d34_r5(X2);
            o[0] = (-24*(pow(m, 6)*p - 6*pow(m, 5)*p + 5*pow(m, 5) + 15*pow(m, 4)*p - 20*pow(m, 3)*p + 15*pow(m, 2)*p - 6*m*p + p)/(pow(m, 5)*pow(m - 1, 6))) + (pow(p, 6))*R0_5 + (-pow(p, 6))*Q1_5;  // mu_mu_mu_mu_mu_mu
            o[1] = (6/pow(m, 4)) + (5*pow(p, 4))*R0_4 + (m*pow(p, 5))*R0_5 + (5*pow(p, 4))*Q1_4 + (-pow(p, 5)*(m - 1))*Q1_5;  // mu_mu_mu_mu_mu_phi
            o[2] = (12*pow(p, 2))*R0_3 + (8*m*pow(p, 3))*R0_4 + (pow(m, 2)*pow(p, 4))*R0_5 + (-12*pow(p, 2))*Q1_3 + (8*pow(p, 3)*(m - 1))*Q1_4 + (-pow(p, 4)*pow(m - 1, 2))*Q1_5;  // mu_mu_mu_mu_phi_phi
            o[3] = (6)*R0_2 + (18*m*p)*R0_3 + (9*pow(m, 2)*pow(p, 2))*R0_4 + (pow(m, 3)*pow(p, 3))*R0_5 + (6)*Q1_2 + (-18*p*(m - 1))*Q1_3 + (9*pow(p, 2)*pow(m - 1, 2))*Q1_4 + (-pow(p, 3)*pow(m - 1, 3))*Q1_5;  // mu_mu_mu_phi_phi_phi
            o[4] = (12*pow(m, 2))*R0_3 + (8*pow(m, 3)*p)*R0_4 + (pow(m, 4)*pow(p, 2))*R0_5 + (-12*pow(m - 1, 2))*Q1_3 + (8*p*pow(m - 1, 3))*Q1_4 + (-pow(p, 2)*pow(m - 1, 4))*Q1_5;  // mu_mu_phi_phi_phi_phi
            o[5] = (12*pow(p, 2))*R0_3 + (8*m*pow(p, 3))*R0_4 + (pow(m, 2)*pow(p, 4))*R0_5 + (-12*pow(p, 2))*Q1_3 + (8*pow(p, 3)*(m - 1))*Q1_4 + (-pow(p, 4)*pow(m - 1, 2))*Q1_5;  // phi_phi_mu_mu_mu_mu
            o[6] = (6)*R0_2 + (18*m*p)*R0_3 + (9*pow(m, 2)*pow(p, 2))*R0_4 + (pow(m, 3)*pow(p, 3))*R0_5 + (6)*Q1_2 + (-18*p*(m - 1))*Q1_3 + (9*pow(p, 2)*pow(m - 1, 2))*Q1_4 + (-pow(p, 3)*pow(m - 1, 3))*Q1_5;  // phi_phi_mu_mu_mu_phi
            o[7] = (12*pow(m, 2))*R0_3 + (8*pow(m, 3)*p)*R0_4 + (pow(m, 4)*pow(p, 2))*R0_5 + (-12*pow(m - 1, 2))*Q1_3 + (8*p*pow(m - 1, 3))*Q1_4 + (-pow(p, 2)*pow(m - 1, 4))*Q1_5;  // phi_phi_mu_mu_phi_phi
            o[8] = (6/pow(p, 4)) + (5*pow(m, 4))*R0_4 + (pow(m, 5)*p)*R0_5 + (5*pow(m - 1, 4))*Q1_4 + (-p*pow(m - 1, 5))*Q1_5;  // phi_phi_mu_phi_phi_phi
            o[9] = (-24*(m*p - p + 5)/pow(p, 6)) + (pow(m, 6))*R0_5 + (-pow(m - 1, 6))*Q1_5 + (-1)*R2_5;  // phi_phi_phi_phi_phi_phi
            o[10] = (-24/(p*pow(m - 1, 5))) + (4*pow(p, 4))*R0_4 + (m*pow(p, 5))*R0_5 + (4*pow(p, 4))*Q1_4 + (-pow(p, 5)*(m - 1))*Q1_5;  // mu_phi_mu_mu_mu_mu
            o[11] = (-6/(pow(p, 2)*pow(m - 1, 4))) + (9*pow(p, 2))*R0_3 + (7*m*pow(p, 3))*R0_4 + (pow(m, 2)*pow(p, 4))*R0_5 + (-9*pow(p, 2))*Q1_3 + (7*pow(p, 3)*(m - 1))*Q1_4 + (-pow(p, 4)*pow(m - 1, 2))*Q1_5;  // mu_phi_mu_mu_mu_phi
            o[12] = (-4/(pow(p, 3)*pow(m - 1, 3))) + (4)*R0_2 + (14*m*p)*R0_3 + (8*pow(m, 2)*pow(p, 2))*R0_4 + (pow(m, 3)*pow(p, 3))*R0_5 + (4)*Q1_2 + (-14*p*(m - 1))*Q1_3 + (8*pow(p, 2)*pow(m - 1, 2))*Q1_4 + (-pow(p, 3)*pow(m - 1, 3))*Q1_5;  // mu_phi_mu_mu_phi_phi
            o[13] = (-6/(pow(p, 4)*pow(m - 1, 2))) + (9*pow(m, 2))*R0_3 + (7*pow(m, 3)*p)*R0_4 + (pow(m, 4)*pow(p, 2))*R0_5 + (-9*pow(m - 1, 2))*Q1_3 + (7*p*pow(m - 1, 3))*Q1_4 + (-pow(p, 2)*pow(m - 1, 4))*Q1_5;  // mu_phi_mu_phi_phi_phi
            o[14] = (-24/(pow(p, 5)*(m - 1))) + (4*pow(m, 4))*R0_4 + (pow(m, 5)*p)*R0_5 + (4*pow(m - 1, 4))*Q1_4 + (-p*pow(m - 1, 5))*Q1_5;  // mu_phi_phi_phi_phi_phi
            break;
        }
        case 6: {
            const double Q1_2 = R::psigamma(1.0 + X1, 2);
            const double Q1_3 = R::psigamma(1.0 + X1, 3);
            const double Q1_4 = R::psigamma(1.0 + X1, 4);
            const double Q1_5 = R::psigamma(1.0 + X1, 5);
            const double Q2_5 = R::psigamma(1.0 + X2, 5);
            const double R0_2 = d7::psi2_rest(X0);
            const double R0_3 = d7::psi3_rest(X0);
            const double R0_4 = d34_r4(X0);
            const double R0_5 = d34_r5(X0);
            o[0] = (-24*(pow(m, 6)*p - 6*pow(m, 5)*p + 5*pow(m, 5) + 15*pow(m, 4)*p - 20*pow(m, 3)*p + 15*pow(m, 2)*p - 6*m*p + p)/(pow(m, 5)*pow(m - 1, 6))) + (pow(p, 6))*R0_5 + (-pow(p, 6))*Q1_5;  // mu_mu_mu_mu_mu_mu
            o[1] = (6/pow(m, 4)) + (5*pow(p, 4))*R0_4 + (m*pow(p, 5))*R0_5 + (5*pow(p, 4))*Q1_4 + (-pow(p, 5)*(m - 1))*Q1_5;  // mu_mu_mu_mu_mu_phi
            o[2] = (12*pow(p, 2))*R0_3 + (8*m*pow(p, 3))*R0_4 + (pow(m, 2)*pow(p, 4))*R0_5 + (-12*pow(p, 2))*Q1_3 + (8*pow(p, 3)*(m - 1))*Q1_4 + (-pow(p, 4)*pow(m - 1, 2))*Q1_5;  // mu_mu_mu_mu_phi_phi
            o[3] = (6)*R0_2 + (18*m*p)*R0_3 + (9*pow(m, 2)*pow(p, 2))*R0_4 + (pow(m, 3)*pow(p, 3))*R0_5 + (6)*Q1_2 + (-18*p*(m - 1))*Q1_3 + (9*pow(p, 2)*pow(m - 1, 2))*Q1_4 + (-pow(p, 3)*pow(m - 1, 3))*Q1_5;  // mu_mu_mu_phi_phi_phi
            o[4] = (12*pow(m, 2))*R0_3 + (8*pow(m, 3)*p)*R0_4 + (pow(m, 4)*pow(p, 2))*R0_5 + (-12*pow(m - 1, 2))*Q1_3 + (8*p*pow(m - 1, 3))*Q1_4 + (-pow(p, 2)*pow(m - 1, 4))*Q1_5;  // mu_mu_phi_phi_phi_phi
            o[5] = (12*pow(p, 2))*R0_3 + (8*m*pow(p, 3))*R0_4 + (pow(m, 2)*pow(p, 4))*R0_5 + (-12*pow(p, 2))*Q1_3 + (8*pow(p, 3)*(m - 1))*Q1_4 + (-pow(p, 4)*pow(m - 1, 2))*Q1_5;  // phi_phi_mu_mu_mu_mu
            o[6] = (6)*R0_2 + (18*m*p)*R0_3 + (9*pow(m, 2)*pow(p, 2))*R0_4 + (pow(m, 3)*pow(p, 3))*R0_5 + (6)*Q1_2 + (-18*p*(m - 1))*Q1_3 + (9*pow(p, 2)*pow(m - 1, 2))*Q1_4 + (-pow(p, 3)*pow(m - 1, 3))*Q1_5;  // phi_phi_mu_mu_mu_phi
            o[7] = (12*pow(m, 2))*R0_3 + (8*pow(m, 3)*p)*R0_4 + (pow(m, 4)*pow(p, 2))*R0_5 + (-12*pow(m - 1, 2))*Q1_3 + (8*p*pow(m - 1, 3))*Q1_4 + (-pow(p, 2)*pow(m - 1, 4))*Q1_5;  // phi_phi_mu_mu_phi_phi
            o[8] = (6/pow(p, 4)) + (5*pow(m, 4))*R0_4 + (pow(m, 5)*p)*R0_5 + (5*pow(m - 1, 4))*Q1_4 + (-p*pow(m - 1, 5))*Q1_5;  // phi_phi_mu_phi_phi_phi
            o[9] = (-24*m/pow(p, 5)) + (pow(m, 6))*R0_5 + (-pow(m - 1, 6))*Q1_5 + (1)*Q2_5;  // phi_phi_phi_phi_phi_phi
            o[10] = (-24/(p*pow(m - 1, 5))) + (4*pow(p, 4))*R0_4 + (m*pow(p, 5))*R0_5 + (4*pow(p, 4))*Q1_4 + (-pow(p, 5)*(m - 1))*Q1_5;  // mu_phi_mu_mu_mu_mu
            o[11] = (-6/(pow(p, 2)*pow(m - 1, 4))) + (9*pow(p, 2))*R0_3 + (7*m*pow(p, 3))*R0_4 + (pow(m, 2)*pow(p, 4))*R0_5 + (-9*pow(p, 2))*Q1_3 + (7*pow(p, 3)*(m - 1))*Q1_4 + (-pow(p, 4)*pow(m - 1, 2))*Q1_5;  // mu_phi_mu_mu_mu_phi
            o[12] = (-4/(pow(p, 3)*pow(m - 1, 3))) + (4)*R0_2 + (14*m*p)*R0_3 + (8*pow(m, 2)*pow(p, 2))*R0_4 + (pow(m, 3)*pow(p, 3))*R0_5 + (4)*Q1_2 + (-14*p*(m - 1))*Q1_3 + (8*pow(p, 2)*pow(m - 1, 2))*Q1_4 + (-pow(p, 3)*pow(m - 1, 3))*Q1_5;  // mu_phi_mu_mu_phi_phi
            o[13] = (-6/(pow(p, 4)*pow(m - 1, 2))) + (9*pow(m, 2))*R0_3 + (7*pow(m, 3)*p)*R0_4 + (pow(m, 4)*pow(p, 2))*R0_5 + (-9*pow(m - 1, 2))*Q1_3 + (7*p*pow(m - 1, 3))*Q1_4 + (-pow(p, 2)*pow(m - 1, 4))*Q1_5;  // mu_phi_mu_phi_phi_phi
            o[14] = (-24/(pow(p, 5)*(m - 1))) + (4*pow(m, 4))*R0_4 + (pow(m, 5)*p)*R0_5 + (4*pow(m - 1, 4))*Q1_4 + (-p*pow(m - 1, 5))*Q1_5;  // mu_phi_phi_phi_phi_phi
            break;
        }
        case 1: {
            const double Q0_2 = R::psigamma(1.0 + X0, 2);
            const double Q0_3 = R::psigamma(1.0 + X0, 3);
            const double Q0_4 = R::psigamma(1.0 + X0, 4);
            const double Q0_5 = R::psigamma(1.0 + X0, 5);
            const double R1_2 = d7::psi2_rest(X1);
            const double R1_3 = d7::psi3_rest(X1);
            const double R1_4 = d34_r4(X1);
            const double R1_5 = d34_r5(X1);
            const double R2_5 = d34_r5(X2);
            o[0] = (24*(pow(m, 6)*p - 5*pow(m, 5) + 25*pow(m, 4) - 50*pow(m, 3) + 50*pow(m, 2) - 25*m + 5)/(pow(m, 6)*pow(m - 1, 5))) + (-pow(p, 6))*Q0_5 + (pow(p, 6))*R1_5;  // mu_mu_mu_mu_mu_mu
            o[1] = (-6/pow(m - 1, 4)) + (-5*pow(p, 4))*Q0_4 + (-m*pow(p, 5))*Q0_5 + (-5*pow(p, 4))*R1_4 + (pow(p, 5)*(m - 1))*R1_5;  // mu_mu_mu_mu_mu_phi
            o[2] = (-12*pow(p, 2))*Q0_3 + (-8*m*pow(p, 3))*Q0_4 + (-pow(m, 2)*pow(p, 4))*Q0_5 + (12*pow(p, 2))*R1_3 + (-8*pow(p, 3)*(m - 1))*R1_4 + (pow(p, 4)*pow(m - 1, 2))*R1_5;  // mu_mu_mu_mu_phi_phi
            o[3] = (-6)*Q0_2 + (-18*m*p)*Q0_3 + (-9*pow(m, 2)*pow(p, 2))*Q0_4 + (-pow(m, 3)*pow(p, 3))*Q0_5 + (-6)*R1_2 + (18*p*(m - 1))*R1_3 + (-9*pow(p, 2)*pow(m - 1, 2))*R1_4 + (pow(p, 3)*pow(m - 1, 3))*R1_5;  // mu_mu_mu_phi_phi_phi
            o[4] = (-12*pow(m, 2))*Q0_3 + (-8*pow(m, 3)*p)*Q0_4 + (-pow(m, 4)*pow(p, 2))*Q0_5 + (12*pow(m - 1, 2))*R1_3 + (-8*p*pow(m - 1, 3))*R1_4 + (pow(p, 2)*pow(m - 1, 4))*R1_5;  // mu_mu_phi_phi_phi_phi
            o[5] = (-12*pow(p, 2))*Q0_3 + (-8*m*pow(p, 3))*Q0_4 + (-pow(m, 2)*pow(p, 4))*Q0_5 + (12*pow(p, 2))*R1_3 + (-8*pow(p, 3)*(m - 1))*R1_4 + (pow(p, 4)*pow(m - 1, 2))*R1_5;  // phi_phi_mu_mu_mu_mu
            o[6] = (-6)*Q0_2 + (-18*m*p)*Q0_3 + (-9*pow(m, 2)*pow(p, 2))*Q0_4 + (-pow(m, 3)*pow(p, 3))*Q0_5 + (-6)*R1_2 + (18*p*(m - 1))*R1_3 + (-9*pow(p, 2)*pow(m - 1, 2))*R1_4 + (pow(p, 3)*pow(m - 1, 3))*R1_5;  // phi_phi_mu_mu_mu_phi
            o[7] = (-12*pow(m, 2))*Q0_3 + (-8*pow(m, 3)*p)*Q0_4 + (-pow(m, 4)*pow(p, 2))*Q0_5 + (12*pow(m - 1, 2))*R1_3 + (-8*p*pow(m - 1, 3))*R1_4 + (pow(p, 2)*pow(m - 1, 4))*R1_5;  // phi_phi_mu_mu_phi_phi
            o[8] = (-6/pow(p, 4)) + (-5*pow(m, 4))*Q0_4 + (-pow(m, 5)*p)*Q0_5 + (-5*pow(m - 1, 4))*R1_4 + (p*pow(m - 1, 5))*R1_5;  // phi_phi_mu_phi_phi_phi
            o[9] = (24*(m*p - 5)/pow(p, 6)) + (-pow(m, 6))*Q0_5 + (pow(m - 1, 6))*R1_5 + (-1)*R2_5;  // phi_phi_phi_phi_phi_phi
            o[10] = (-24/(pow(m, 5)*p)) + (-4*pow(p, 4))*Q0_4 + (-m*pow(p, 5))*Q0_5 + (-4*pow(p, 4))*R1_4 + (pow(p, 5)*(m - 1))*R1_5;  // mu_phi_mu_mu_mu_mu
            o[11] = (-6/(pow(m, 4)*pow(p, 2))) + (-9*pow(p, 2))*Q0_3 + (-7*m*pow(p, 3))*Q0_4 + (-pow(m, 2)*pow(p, 4))*Q0_5 + (9*pow(p, 2))*R1_3 + (-7*pow(p, 3)*(m - 1))*R1_4 + (pow(p, 4)*pow(m - 1, 2))*R1_5;  // mu_phi_mu_mu_mu_phi
            o[12] = (-4/(pow(m, 3)*pow(p, 3))) + (-4)*Q0_2 + (-14*m*p)*Q0_3 + (-8*pow(m, 2)*pow(p, 2))*Q0_4 + (-pow(m, 3)*pow(p, 3))*Q0_5 + (-4)*R1_2 + (14*p*(m - 1))*R1_3 + (-8*pow(p, 2)*pow(m - 1, 2))*R1_4 + (pow(p, 3)*pow(m - 1, 3))*R1_5;  // mu_phi_mu_mu_phi_phi
            o[13] = (-6/(pow(m, 2)*pow(p, 4))) + (-9*pow(m, 2))*Q0_3 + (-7*pow(m, 3)*p)*Q0_4 + (-pow(m, 4)*pow(p, 2))*Q0_5 + (9*pow(m - 1, 2))*R1_3 + (-7*p*pow(m - 1, 3))*R1_4 + (pow(p, 2)*pow(m - 1, 4))*R1_5;  // mu_phi_mu_phi_phi_phi
            o[14] = (-24/(m*pow(p, 5))) + (-4*pow(m, 4))*Q0_4 + (-pow(m, 5)*p)*Q0_5 + (-4*pow(m - 1, 4))*R1_4 + (p*pow(m - 1, 5))*R1_5;  // mu_phi_phi_phi_phi_phi
            break;
        }
        case 5: {
            const double Q0_2 = R::psigamma(1.0 + X0, 2);
            const double Q0_3 = R::psigamma(1.0 + X0, 3);
            const double Q0_4 = R::psigamma(1.0 + X0, 4);
            const double Q0_5 = R::psigamma(1.0 + X0, 5);
            const double Q2_5 = R::psigamma(1.0 + X2, 5);
            const double R1_2 = d7::psi2_rest(X1);
            const double R1_3 = d7::psi3_rest(X1);
            const double R1_4 = d34_r4(X1);
            const double R1_5 = d34_r5(X1);
            o[0] = (24*(pow(m, 6)*p - 5*pow(m, 5) + 25*pow(m, 4) - 50*pow(m, 3) + 50*pow(m, 2) - 25*m + 5)/(pow(m, 6)*pow(m - 1, 5))) + (-pow(p, 6))*Q0_5 + (pow(p, 6))*R1_5;  // mu_mu_mu_mu_mu_mu
            o[1] = (-6/pow(m - 1, 4)) + (-5*pow(p, 4))*Q0_4 + (-m*pow(p, 5))*Q0_5 + (-5*pow(p, 4))*R1_4 + (pow(p, 5)*(m - 1))*R1_5;  // mu_mu_mu_mu_mu_phi
            o[2] = (-12*pow(p, 2))*Q0_3 + (-8*m*pow(p, 3))*Q0_4 + (-pow(m, 2)*pow(p, 4))*Q0_5 + (12*pow(p, 2))*R1_3 + (-8*pow(p, 3)*(m - 1))*R1_4 + (pow(p, 4)*pow(m - 1, 2))*R1_5;  // mu_mu_mu_mu_phi_phi
            o[3] = (-6)*Q0_2 + (-18*m*p)*Q0_3 + (-9*pow(m, 2)*pow(p, 2))*Q0_4 + (-pow(m, 3)*pow(p, 3))*Q0_5 + (-6)*R1_2 + (18*p*(m - 1))*R1_3 + (-9*pow(p, 2)*pow(m - 1, 2))*R1_4 + (pow(p, 3)*pow(m - 1, 3))*R1_5;  // mu_mu_mu_phi_phi_phi
            o[4] = (-12*pow(m, 2))*Q0_3 + (-8*pow(m, 3)*p)*Q0_4 + (-pow(m, 4)*pow(p, 2))*Q0_5 + (12*pow(m - 1, 2))*R1_3 + (-8*p*pow(m - 1, 3))*R1_4 + (pow(p, 2)*pow(m - 1, 4))*R1_5;  // mu_mu_phi_phi_phi_phi
            o[5] = (-12*pow(p, 2))*Q0_3 + (-8*m*pow(p, 3))*Q0_4 + (-pow(m, 2)*pow(p, 4))*Q0_5 + (12*pow(p, 2))*R1_3 + (-8*pow(p, 3)*(m - 1))*R1_4 + (pow(p, 4)*pow(m - 1, 2))*R1_5;  // phi_phi_mu_mu_mu_mu
            o[6] = (-6)*Q0_2 + (-18*m*p)*Q0_3 + (-9*pow(m, 2)*pow(p, 2))*Q0_4 + (-pow(m, 3)*pow(p, 3))*Q0_5 + (-6)*R1_2 + (18*p*(m - 1))*R1_3 + (-9*pow(p, 2)*pow(m - 1, 2))*R1_4 + (pow(p, 3)*pow(m - 1, 3))*R1_5;  // phi_phi_mu_mu_mu_phi
            o[7] = (-12*pow(m, 2))*Q0_3 + (-8*pow(m, 3)*p)*Q0_4 + (-pow(m, 4)*pow(p, 2))*Q0_5 + (12*pow(m - 1, 2))*R1_3 + (-8*p*pow(m - 1, 3))*R1_4 + (pow(p, 2)*pow(m - 1, 4))*R1_5;  // phi_phi_mu_mu_phi_phi
            o[8] = (-6/pow(p, 4)) + (-5*pow(m, 4))*Q0_4 + (-pow(m, 5)*p)*Q0_5 + (-5*pow(m - 1, 4))*R1_4 + (p*pow(m - 1, 5))*R1_5;  // phi_phi_mu_phi_phi_phi
            o[9] = (24*(m - 1)/pow(p, 5)) + (-pow(m, 6))*Q0_5 + (pow(m - 1, 6))*R1_5 + (1)*Q2_5;  // phi_phi_phi_phi_phi_phi
            o[10] = (-24/(pow(m, 5)*p)) + (-4*pow(p, 4))*Q0_4 + (-m*pow(p, 5))*Q0_5 + (-4*pow(p, 4))*R1_4 + (pow(p, 5)*(m - 1))*R1_5;  // mu_phi_mu_mu_mu_mu
            o[11] = (-6/(pow(m, 4)*pow(p, 2))) + (-9*pow(p, 2))*Q0_3 + (-7*m*pow(p, 3))*Q0_4 + (-pow(m, 2)*pow(p, 4))*Q0_5 + (9*pow(p, 2))*R1_3 + (-7*pow(p, 3)*(m - 1))*R1_4 + (pow(p, 4)*pow(m - 1, 2))*R1_5;  // mu_phi_mu_mu_mu_phi
            o[12] = (-4/(pow(m, 3)*pow(p, 3))) + (-4)*Q0_2 + (-14*m*p)*Q0_3 + (-8*pow(m, 2)*pow(p, 2))*Q0_4 + (-pow(m, 3)*pow(p, 3))*Q0_5 + (-4)*R1_2 + (14*p*(m - 1))*R1_3 + (-8*pow(p, 2)*pow(m - 1, 2))*R1_4 + (pow(p, 3)*pow(m - 1, 3))*R1_5;  // mu_phi_mu_mu_phi_phi
            o[13] = (-6/(pow(m, 2)*pow(p, 4))) + (-9*pow(m, 2))*Q0_3 + (-7*pow(m, 3)*p)*Q0_4 + (-pow(m, 4)*pow(p, 2))*Q0_5 + (9*pow(m - 1, 2))*R1_3 + (-7*p*pow(m - 1, 3))*R1_4 + (pow(p, 2)*pow(m - 1, 4))*R1_5;  // mu_phi_mu_phi_phi_phi
            o[14] = (-24/(m*pow(p, 5))) + (-4*pow(m, 4))*Q0_4 + (-pow(m, 5)*p)*Q0_5 + (-4*pow(m - 1, 4))*R1_4 + (p*pow(m - 1, 5))*R1_5;  // mu_phi_phi_phi_phi_phi
            break;
        }
        case 3: {
            const double Q0_2 = R::psigamma(1.0 + X0, 2);
            const double Q0_3 = R::psigamma(1.0 + X0, 3);
            const double Q0_4 = R::psigamma(1.0 + X0, 4);
            const double Q0_5 = R::psigamma(1.0 + X0, 5);
            const double Q1_2 = R::psigamma(1.0 + X1, 2);
            const double Q1_3 = R::psigamma(1.0 + X1, 3);
            const double Q1_4 = R::psigamma(1.0 + X1, 4);
            const double Q1_5 = R::psigamma(1.0 + X1, 5);
            const double R2_5 = d34_r5(X2);
            o[0] = (-120*(2*pow(m, 2) - 2*m + 1)*(pow(m, 4) - 2*pow(m, 3) + 5*pow(m, 2) - 4*m + 1)/(pow(m, 6)*pow(m - 1, 6))) + (-pow(p, 6))*Q0_5 + (-pow(p, 6))*Q1_5;  // mu_mu_mu_mu_mu_mu
            o[1] = (-5*pow(p, 4))*Q0_4 + (-m*pow(p, 5))*Q0_5 + (5*pow(p, 4))*Q1_4 + (-pow(p, 5)*(m - 1))*Q1_5;  // mu_mu_mu_mu_mu_phi
            o[2] = (-12*pow(p, 2))*Q0_3 + (-8*m*pow(p, 3))*Q0_4 + (-pow(m, 2)*pow(p, 4))*Q0_5 + (-12*pow(p, 2))*Q1_3 + (8*pow(p, 3)*(m - 1))*Q1_4 + (-pow(p, 4)*pow(m - 1, 2))*Q1_5;  // mu_mu_mu_mu_phi_phi
            o[3] = (-6)*Q0_2 + (-18*m*p)*Q0_3 + (-9*pow(m, 2)*pow(p, 2))*Q0_4 + (-pow(m, 3)*pow(p, 3))*Q0_5 + (6)*Q1_2 + (-18*p*(m - 1))*Q1_3 + (9*pow(p, 2)*pow(m - 1, 2))*Q1_4 + (-pow(p, 3)*pow(m - 1, 3))*Q1_5;  // mu_mu_mu_phi_phi_phi
            o[4] = (-12*pow(m, 2))*Q0_3 + (-8*pow(m, 3)*p)*Q0_4 + (-pow(m, 4)*pow(p, 2))*Q0_5 + (-12*pow(m - 1, 2))*Q1_3 + (8*p*pow(m - 1, 3))*Q1_4 + (-pow(p, 2)*pow(m - 1, 4))*Q1_5;  // mu_mu_phi_phi_phi_phi
            o[5] = (-12*pow(p, 2))*Q0_3 + (-8*m*pow(p, 3))*Q0_4 + (-pow(m, 2)*pow(p, 4))*Q0_5 + (-12*pow(p, 2))*Q1_3 + (8*pow(p, 3)*(m - 1))*Q1_4 + (-pow(p, 4)*pow(m - 1, 2))*Q1_5;  // phi_phi_mu_mu_mu_mu
            o[6] = (-6)*Q0_2 + (-18*m*p)*Q0_3 + (-9*pow(m, 2)*pow(p, 2))*Q0_4 + (-pow(m, 3)*pow(p, 3))*Q0_5 + (6)*Q1_2 + (-18*p*(m - 1))*Q1_3 + (9*pow(p, 2)*pow(m - 1, 2))*Q1_4 + (-pow(p, 3)*pow(m - 1, 3))*Q1_5;  // phi_phi_mu_mu_mu_phi
            o[7] = (-12*pow(m, 2))*Q0_3 + (-8*pow(m, 3)*p)*Q0_4 + (-pow(m, 4)*pow(p, 2))*Q0_5 + (-12*pow(m - 1, 2))*Q1_3 + (8*p*pow(m - 1, 3))*Q1_4 + (-pow(p, 2)*pow(m - 1, 4))*Q1_5;  // phi_phi_mu_mu_phi_phi
            o[8] = (-5*pow(m, 4))*Q0_4 + (-pow(m, 5)*p)*Q0_5 + (5*pow(m - 1, 4))*Q1_4 + (-p*pow(m - 1, 5))*Q1_5;  // phi_phi_mu_phi_phi_phi
            o[9] = (24*(p - 10)/pow(p, 6)) + (-pow(m, 6))*Q0_5 + (-pow(m - 1, 6))*Q1_5 + (-1)*R2_5;  // phi_phi_phi_phi_phi_phi
            o[10] = (-24*(2*m - 1)*(pow(m, 4) - 2*pow(m, 3) + 4*pow(m, 2) - 3*m + 1)/(pow(m, 5)*p*pow(m - 1, 5))) + (-4*pow(p, 4))*Q0_4 + (-m*pow(p, 5))*Q0_5 + (4*pow(p, 4))*Q1_4 + (-pow(p, 5)*(m - 1))*Q1_5;  // mu_phi_mu_mu_mu_mu
            o[11] = (-6*(2*pow(m, 4) - 4*pow(m, 3) + 6*pow(m, 2) - 4*m + 1)/(pow(m, 4)*pow(p, 2)*pow(m - 1, 4))) + (-9*pow(p, 2))*Q0_3 + (-7*m*pow(p, 3))*Q0_4 + (-pow(m, 2)*pow(p, 4))*Q0_5 + (-9*pow(p, 2))*Q1_3 + (7*pow(p, 3)*(m - 1))*Q1_4 + (-pow(p, 4)*pow(m - 1, 2))*Q1_5;  // mu_phi_mu_mu_mu_phi
            o[12] = (-4*(2*m - 1)*(pow(m, 2) - m + 1)/(pow(m, 3)*pow(p, 3)*pow(m - 1, 3))) + (-4)*Q0_2 + (-14*m*p)*Q0_3 + (-8*pow(m, 2)*pow(p, 2))*Q0_4 + (-pow(m, 3)*pow(p, 3))*Q0_5 + (4)*Q1_2 + (-14*p*(m - 1))*Q1_3 + (8*pow(p, 2)*pow(m - 1, 2))*Q1_4 + (-pow(p, 3)*pow(m - 1, 3))*Q1_5;  // mu_phi_mu_mu_phi_phi
            o[13] = (-6*(2*pow(m, 2) - 2*m + 1)/(pow(m, 2)*pow(p, 4)*pow(m - 1, 2))) + (-9*pow(m, 2))*Q0_3 + (-7*pow(m, 3)*p)*Q0_4 + (-pow(m, 4)*pow(p, 2))*Q0_5 + (-9*pow(m - 1, 2))*Q1_3 + (7*p*pow(m - 1, 3))*Q1_4 + (-pow(p, 2)*pow(m - 1, 4))*Q1_5;  // mu_phi_mu_phi_phi_phi
            o[14] = (-24*(2*m - 1)/(m*pow(p, 5)*(m - 1))) + (-4*pow(m, 4))*Q0_4 + (-pow(m, 5)*p)*Q0_5 + (4*pow(m - 1, 4))*Q1_4 + (-p*pow(m - 1, 5))*Q1_5;  // mu_phi_phi_phi_phi_phi
            break;
        }
        case 7: {
            const double Q0_2 = R::psigamma(1.0 + X0, 2);
            const double Q0_3 = R::psigamma(1.0 + X0, 3);
            const double Q0_4 = R::psigamma(1.0 + X0, 4);
            const double Q0_5 = R::psigamma(1.0 + X0, 5);
            const double Q1_2 = R::psigamma(1.0 + X1, 2);
            const double Q1_3 = R::psigamma(1.0 + X1, 3);
            const double Q1_4 = R::psigamma(1.0 + X1, 4);
            const double Q1_5 = R::psigamma(1.0 + X1, 5);
            const double Q2_5 = R::psigamma(1.0 + X2, 5);
            o[0] = (-120*(2*pow(m, 2) - 2*m + 1)*(pow(m, 4) - 2*pow(m, 3) + 5*pow(m, 2) - 4*m + 1)/(pow(m, 6)*pow(m - 1, 6))) + (-pow(p, 6))*Q0_5 + (-pow(p, 6))*Q1_5;  // mu_mu_mu_mu_mu_mu
            o[1] = (-5*pow(p, 4))*Q0_4 + (-m*pow(p, 5))*Q0_5 + (5*pow(p, 4))*Q1_4 + (-pow(p, 5)*(m - 1))*Q1_5;  // mu_mu_mu_mu_mu_phi
            o[2] = (-12*pow(p, 2))*Q0_3 + (-8*m*pow(p, 3))*Q0_4 + (-pow(m, 2)*pow(p, 4))*Q0_5 + (-12*pow(p, 2))*Q1_3 + (8*pow(p, 3)*(m - 1))*Q1_4 + (-pow(p, 4)*pow(m - 1, 2))*Q1_5;  // mu_mu_mu_mu_phi_phi
            o[3] = (-6)*Q0_2 + (-18*m*p)*Q0_3 + (-9*pow(m, 2)*pow(p, 2))*Q0_4 + (-pow(m, 3)*pow(p, 3))*Q0_5 + (6)*Q1_2 + (-18*p*(m - 1))*Q1_3 + (9*pow(p, 2)*pow(m - 1, 2))*Q1_4 + (-pow(p, 3)*pow(m - 1, 3))*Q1_5;  // mu_mu_mu_phi_phi_phi
            o[4] = (-12*pow(m, 2))*Q0_3 + (-8*pow(m, 3)*p)*Q0_4 + (-pow(m, 4)*pow(p, 2))*Q0_5 + (-12*pow(m - 1, 2))*Q1_3 + (8*p*pow(m - 1, 3))*Q1_4 + (-pow(p, 2)*pow(m - 1, 4))*Q1_5;  // mu_mu_phi_phi_phi_phi
            o[5] = (-12*pow(p, 2))*Q0_3 + (-8*m*pow(p, 3))*Q0_4 + (-pow(m, 2)*pow(p, 4))*Q0_5 + (-12*pow(p, 2))*Q1_3 + (8*pow(p, 3)*(m - 1))*Q1_4 + (-pow(p, 4)*pow(m - 1, 2))*Q1_5;  // phi_phi_mu_mu_mu_mu
            o[6] = (-6)*Q0_2 + (-18*m*p)*Q0_3 + (-9*pow(m, 2)*pow(p, 2))*Q0_4 + (-pow(m, 3)*pow(p, 3))*Q0_5 + (6)*Q1_2 + (-18*p*(m - 1))*Q1_3 + (9*pow(p, 2)*pow(m - 1, 2))*Q1_4 + (-pow(p, 3)*pow(m - 1, 3))*Q1_5;  // phi_phi_mu_mu_mu_phi
            o[7] = (-12*pow(m, 2))*Q0_3 + (-8*pow(m, 3)*p)*Q0_4 + (-pow(m, 4)*pow(p, 2))*Q0_5 + (-12*pow(m - 1, 2))*Q1_3 + (8*p*pow(m - 1, 3))*Q1_4 + (-pow(p, 2)*pow(m - 1, 4))*Q1_5;  // phi_phi_mu_mu_phi_phi
            o[8] = (-5*pow(m, 4))*Q0_4 + (-pow(m, 5)*p)*Q0_5 + (5*pow(m - 1, 4))*Q1_4 + (-p*pow(m - 1, 5))*Q1_5;  // phi_phi_mu_phi_phi_phi
            o[9] = (-120/pow(p, 6)) + (-pow(m, 6))*Q0_5 + (-pow(m - 1, 6))*Q1_5 + (1)*Q2_5;  // phi_phi_phi_phi_phi_phi
            o[10] = (-24*(2*m - 1)*(pow(m, 4) - 2*pow(m, 3) + 4*pow(m, 2) - 3*m + 1)/(pow(m, 5)*p*pow(m - 1, 5))) + (-4*pow(p, 4))*Q0_4 + (-m*pow(p, 5))*Q0_5 + (4*pow(p, 4))*Q1_4 + (-pow(p, 5)*(m - 1))*Q1_5;  // mu_phi_mu_mu_mu_mu
            o[11] = (-6*(2*pow(m, 4) - 4*pow(m, 3) + 6*pow(m, 2) - 4*m + 1)/(pow(m, 4)*pow(p, 2)*pow(m - 1, 4))) + (-9*pow(p, 2))*Q0_3 + (-7*m*pow(p, 3))*Q0_4 + (-pow(m, 2)*pow(p, 4))*Q0_5 + (-9*pow(p, 2))*Q1_3 + (7*pow(p, 3)*(m - 1))*Q1_4 + (-pow(p, 4)*pow(m - 1, 2))*Q1_5;  // mu_phi_mu_mu_mu_phi
            o[12] = (-4*(2*m - 1)*(pow(m, 2) - m + 1)/(pow(m, 3)*pow(p, 3)*pow(m - 1, 3))) + (-4)*Q0_2 + (-14*m*p)*Q0_3 + (-8*pow(m, 2)*pow(p, 2))*Q0_4 + (-pow(m, 3)*pow(p, 3))*Q0_5 + (4)*Q1_2 + (-14*p*(m - 1))*Q1_3 + (8*pow(p, 2)*pow(m - 1, 2))*Q1_4 + (-pow(p, 3)*pow(m - 1, 3))*Q1_5;  // mu_phi_mu_mu_phi_phi
            o[13] = (-6*(2*pow(m, 2) - 2*m + 1)/(pow(m, 2)*pow(p, 4)*pow(m - 1, 2))) + (-9*pow(m, 2))*Q0_3 + (-7*pow(m, 3)*p)*Q0_4 + (-pow(m, 4)*pow(p, 2))*Q0_5 + (-9*pow(m - 1, 2))*Q1_3 + (7*p*pow(m - 1, 3))*Q1_4 + (-pow(p, 2)*pow(m - 1, 4))*Q1_5;  // mu_phi_mu_phi_phi_phi
            o[14] = (-24*(2*m - 1)/(m*pow(p, 5)*(m - 1))) + (-4*pow(m, 4))*Q0_4 + (-pow(m, 5)*p)*Q0_5 + (4*pow(m - 1, 4))*Q1_4 + (-p*pow(m - 1, 5))*Q1_5;  // mu_phi_phi_phi_phi_phi
            break;
        }
        }
    });
}
