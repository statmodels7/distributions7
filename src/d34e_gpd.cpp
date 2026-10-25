// The third and fourth derivatives of the expected information of the gpd
// family on the parameter scale, in closed form. The order-2 kernel of
// distributions7 states E_ab in closed form; the kernels below differentiate
// that closed form, one order per kernel, and write each component in the
// order of d3expected_names() and d4expected_names() for gpd.

#include <Rcpp.h>
#include "d7_par.h"
#include "d34e_run.h"
using namespace Rcpp;

// [[Rcpp::export]]
List gpd_dexpected3_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
    d34e::Par S(sigma);
    d34e::Par X(xi);
    return d34e::run(y.size(), threads, d7::kMinCheap, 12,
                     [&](std::size_t i, d34e::Row o) {
        const double s = S[i], x = X[i];
        if (x <= -0.5) {
            for (int k = 0; k < 12; ++k) o[k] = NA_REAL;
            return;
        }
        o[0] = 24/(pow(s, 5)*(2*x + 1));   // sigma_sigma_sigma_sigma_sigma
        o[1] = 12/(pow(s, 4)*pow(2*x + 1, 2));   // sigma_sigma_sigma_sigma_xi
        o[2] = 16/(pow(s, 3)*pow(2*x + 1, 3));   // sigma_sigma_sigma_xi_xi
        o[3] = 48/(pow(s, 2)*pow(2*x + 1, 4));   // sigma_sigma_xi_xi_xi
        o[4] = 0;   // xi_xi_sigma_sigma_sigma
        o[5] = 0;   // xi_xi_sigma_sigma_xi
        o[6] = 0;   // xi_xi_sigma_xi_xi
        o[7] = 12*(4*x + 3)*(8*pow(x, 2) + 12*x + 5)/(pow(x + 1, 4)*pow(2*x + 1, 4));   // xi_xi_xi_xi_xi
        o[8] = 6/(pow(s, 4)*(x + 1)*(2*x + 1));   // sigma_xi_sigma_sigma_sigma
        o[9] = 2*(4*x + 3)/(pow(s, 3)*pow(x + 1, 2)*pow(2*x + 1, 2));   // sigma_xi_sigma_sigma_xi
        o[10] = 2*(12*pow(x, 2) + 18*x + 7)/(pow(s, 2)*pow(x + 1, 3)*pow(2*x + 1, 3));   // sigma_xi_sigma_xi_xi
        o[11] = 6*(4*x + 3)*(8*pow(x, 2) + 12*x + 5)/(s*pow(x + 1, 4)*pow(2*x + 1, 4));   // sigma_xi_xi_xi_xi
    });
}

// [[Rcpp::export]]
List gpd_dexpected4_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
    d34e::Par S(sigma);
    d34e::Par X(xi);
    return d34e::run(y.size(), threads, d7::kMinCheap, 15,
                     [&](std::size_t i, d34e::Row o) {
        const double s = S[i], x = X[i];
        if (x <= -0.5) {
            for (int k = 0; k < 15; ++k) o[k] = NA_REAL;
            return;
        }
        o[0] = -120/(pow(s, 6)*(2*x + 1));   // sigma_sigma_sigma_sigma_sigma_sigma
        o[1] = -48/(pow(s, 5)*pow(2*x + 1, 2));   // sigma_sigma_sigma_sigma_sigma_xi
        o[2] = -48/(pow(s, 4)*pow(2*x + 1, 3));   // sigma_sigma_sigma_sigma_xi_xi
        o[3] = -96/(pow(s, 3)*pow(2*x + 1, 4));   // sigma_sigma_sigma_xi_xi_xi
        o[4] = -384/(pow(s, 2)*pow(2*x + 1, 5));   // sigma_sigma_xi_xi_xi_xi
        o[5] = 0;   // xi_xi_sigma_sigma_sigma_sigma
        o[6] = 0;   // xi_xi_sigma_sigma_sigma_xi
        o[7] = 0;   // xi_xi_sigma_sigma_xi_xi
        o[8] = 0;   // xi_xi_sigma_xi_xi_xi
        o[9] = -48*(80*pow(x, 4) + 240*pow(x, 3) + 280*pow(x, 2) + 150*x + 31)/(pow(x + 1, 5)*pow(2*x + 1, 5));   // xi_xi_xi_xi_xi_xi
        o[10] = -24/(pow(s, 5)*(x + 1)*(2*x + 1));   // sigma_xi_sigma_sigma_sigma_sigma
        o[11] = -6*(4*x + 3)/(pow(s, 4)*pow(x + 1, 2)*pow(2*x + 1, 2));   // sigma_xi_sigma_sigma_sigma_xi
        o[12] = -4*(12*pow(x, 2) + 18*x + 7)/(pow(s, 3)*pow(x + 1, 3)*pow(2*x + 1, 3));   // sigma_xi_sigma_sigma_xi_xi
        o[13] = -6*(4*x + 3)*(8*pow(x, 2) + 12*x + 5)/(pow(s, 2)*pow(x + 1, 4)*pow(2*x + 1, 4));   // sigma_xi_sigma_xi_xi_xi
        o[14] = -24*(80*pow(x, 4) + 240*pow(x, 3) + 280*pow(x, 2) + 150*x + 31)/(s*pow(x + 1, 5)*pow(2*x + 1, 5));   // sigma_xi_xi_xi_xi_xi
    });
}
