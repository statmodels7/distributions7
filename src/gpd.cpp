#include <Rcpp.h>
using namespace Rcpp;

// Generalised Pareto in (sigma, xi):
//   f(y) = (1/sigma) (1 + xi y/sigma)^(-1/xi - 1),
//   l    = -log(sigma) - log(t)/xi - log(t),   t = 1 + xi z,  z = y/sigma.
// Writing w = log(t)/xi removes the apparent singularity at xi = 0, where
// w -> z and the family becomes the exponential.
//
// With u = z/t, the two facts that keep every derivative short are
// t - xi z = 1, so du/dsigma = -z/(sigma t^2), and dt/dxi = z, so
// du/dxi = -z^2/t^2.

// log(1 + xi z)/xi, stable through zero.
static inline double gpd_w(double xi, double z) {
    if (std::fabs(xi) < 1e-8) {
        return z - 0.5 * xi * z * z + xi * xi * z * z * z / 3.0;
    }
    return std::log1p(xi * z) / xi;
}

// dl/dxi, whose limit at xi = 0 is z^2/2 - z.
static inline double gpd_dl_dxi(double xi, double z, double t) {
    if (std::fabs(xi) < 1e-6) {
        double z2 = z * z, z3 = z2 * z;
        return 0.5 * z2 - z + xi * (z3 * (-2.0 / 3.0) + z2);
    }
    return std::log1p(xi * z) / (xi * xi) - (1.0 + 1.0 / xi) * z / t;
}

// [[Rcpp::export]]
NumericVector gpd_logpdf_cpp(NumericVector y, NumericVector sigma,
                             NumericVector xi) {
    int n = y.size();
    NumericVector out(n);
    bool s_s = (sigma.size() == 1), x_s = (xi.size() == 1);
    for (int i = 0; i < n; i++) {
        double s = s_s ? sigma[0] : sigma[i];
        double x = x_s ? xi[0] : xi[i];
        double z = y[i] / s, t = 1.0 + x * z;
        if (y[i] < 0 || t <= 0) { out[i] = R_NegInf; continue; }
        out[i] = -std::log(s) - gpd_w(x, z) - std::log(t);
    }
    return out;
}

// [[Rcpp::export]]
List gpd_gradient_cpp(NumericVector y, NumericVector sigma, NumericVector xi) {
    int n = y.size();
    NumericVector g_s(n), g_x(n);
    bool s_s = (sigma.size() == 1), x_s = (xi.size() == 1);
    for (int i = 0; i < n; i++) {
        double s = s_s ? sigma[0] : sigma[i];
        double x = x_s ? xi[0] : xi[i];
        double z = y[i] / s, t = 1.0 + x * z, u = z / t;
        g_s[i] = ((x + 1.0) * u - 1.0) / s;
        g_x[i] = gpd_dl_dxi(x, z, t);
    }
    return List::create(Named("sigma") = g_s, Named("xi") = g_x);
}

// [[Rcpp::export]]
List gpd_hessian_cpp(NumericVector y, NumericVector sigma, NumericVector xi) {
    int n = y.size();
    NumericVector h_ss(n), h_sx(n), h_xx(n);
    bool s_s = (sigma.size() == 1), x_s = (xi.size() == 1);
    for (int i = 0; i < n; i++) {
        double s = s_s ? sigma[0] : sigma[i];
        double x = x_s ? xi[0] : xi[i];
        double z = y[i] / s, t = 1.0 + x * z, u = z / t;
        double t2 = t * t, z2 = z * z;

        h_ss[i] = (-(x + 1.0) * z / t2 - (x + 1.0) * u + 1.0) / (s * s);
        h_sx[i] = (u - (x + 1.0) * z2 / t2) / s;

        if (std::fabs(x) < 1e-4) {
            // the xi-xi component through its series: the two singular terms
            // of the closed form cancel and what remains is O(1)
            double z3 = z2 * z, z4 = z3 * z;
            h_xx[i] = -2.0 * z3 / 3.0 + z2 + x * (1.5 * z4 - 2.0 * z3);
        } else {
            h_xx[i] = -2.0 * std::log1p(x * z) / (x * x * x)
                    + 2.0 * u / (x * x)
                    + (1.0 + 1.0 / x) * z2 / t2;
        }
    }
    return List::create(Named("sigma_sigma") = h_ss,
                        Named("sigma_xi") = h_sx,
                        Named("xi_xi") = h_xx);
}

// The expected information, closed form for xi > -1/2 (Smith, 1985):
//   I = 1/(1 + 2 xi) * [ 1/sigma^2            1/(sigma (1 + xi))
//                        1/(sigma (1 + xi))   2/(1 + xi)        ]
// Below -1/2 it does not exist and the classical asymptotics do not hold,
// which is a property of the family rather than of this implementation.
// [[Rcpp::export]]
List gpd_expected_hessian_cpp(NumericVector y, NumericVector sigma,
                              NumericVector xi) {
    int n = y.size();
    NumericVector h_ss(n), h_sx(n), h_xx(n);
    bool s_s = (sigma.size() == 1), x_s = (xi.size() == 1);
    for (int i = 0; i < n; i++) {
        double s = s_s ? sigma[0] : sigma[i];
        double x = x_s ? xi[0] : xi[i];
        if (x <= -0.5) {
            h_ss[i] = NA_REAL; h_sx[i] = NA_REAL; h_xx[i] = NA_REAL;
            continue;
        }
        double d = 1.0 + 2.0 * x, om = 1.0 + x;
        h_ss[i] = -1.0 / (d * s * s);
        h_sx[i] = -1.0 / (d * s * om);
        h_xx[i] = -2.0 / (d * om);
    }
    return List::create(Named("sigma_sigma") = h_ss,
                        Named("sigma_xi") = h_sx,
                        Named("xi_xi") = h_xx);
}
