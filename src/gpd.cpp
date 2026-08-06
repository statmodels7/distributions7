#include <Rcpp.h>
using namespace Rcpp;

// Generalized Pareto in (sigma, xi):
//   f(y) = (1/sigma) (1 + xi y/sigma)^(-1/xi - 1),
//   l    = -log(sigma) - log(t)/xi - log(t),   t = 1 + xi z,  z = y/sigma.
// Writing w = log(t)/xi removes the apparent singularity at xi = 0, where
// w -> z and the family becomes the exponential.
//
// With u = z/t, the two facts that keep every derivative short are
// t - xi z = 1, so du/dsigma = -z/(sigma t^2), and dt/dxi = z, so
// du/dxi = -z^2/t^2.

// The xi direction runs through W = log1p(xi z)/xi and its xi-derivatives.
// W itself never cancels: written as log1p(u) * (z/u) with u = xi z, the
// relative error is machine precision for every u > -1 (a guard on |xi|
// alone, with a short series in xi, was wrong here too -- at xi = 1e-8 and
// z = 1e7 the truncation is (xi z)^3/4, four decimals). The derivatives DO
// cancel: at order b the Leibniz terms are of size z xi^-b against an
// answer of size z^(b+1), so the relative cancellation is (xi z)^-b --
// measured, the order-two direct form keeps a floor of 1.8e-11 already at
// xi z = 5e-4, which is what the old |xi| < 1e-4 guard missed. Below
// |xi z| = 0.2 both derivatives are taken from the series of W
// differentiated term by term; forty terms leave (0.2)^39 ~ 2e-28 and the
// direct forms still carry ~5e-15 at the cut, so the two routes agree in
// the overlap.

static inline double gpd_w(double xi, double z) {
    double u = xi * z;
    if (u == 0.0) return z;
    return std::log1p(u) * (z / u);
}

// d/dxi of W: series sum_{k>=1} (-1)^k k (xi z)^(k-1) z^2 / (k+1).
static inline double gpd_w_dxi(double xi, double z, double t) {
    double u = xi * z;
    if (std::fabs(u) < 0.2) {
        double s = 0.0, uk = 1.0, sgn = -1.0;
        for (int k = 1; k <= 40; ++k) {
            s += sgn * k * uk / (k + 1.0);
            uk *= u; sgn = -sgn;
        }
        return s * z * z;
    }
    return (z / t) / xi - std::log1p(u) / (xi * xi);
}

// d^2/dxi^2 of W: series sum_{k>=2} (-1)^k k(k-1) (xi z)^(k-2) z^3 / (k+1).
static inline double gpd_w_d2xi(double xi, double z, double t) {
    double u = xi * z;
    if (std::fabs(u) < 0.2) {
        double s = 0.0, uk = 1.0, sgn = 1.0;
        for (int k = 2; k <= 41; ++k) {
            s += sgn * k * (k - 1.0) * uk / (k + 1.0);
            uk *= u; sgn = -sgn;
        }
        return s * z * z * z;
    }
    return -(z * z / (t * t)) / xi - 2.0 * (z / t) / (xi * xi)
         + 2.0 * std::log1p(u) / (xi * xi * xi);
}

// dl/dxi = -z/t - dW/dxi, whose limit at xi = 0 is z^2/2 - z.
static inline double gpd_dl_dxi(double xi, double z, double t) {
    return -z / t - gpd_w_dxi(xi, z, t);
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

        // l = -log(sigma) - log(t) - W, so l_xixi = z^2/t^2 - d^2 W/dxi^2,
        // the second piece switching to its series where the direct form's
        // terms of size z^2 xi^-2 cancel (see gpd_w_d2xi above)
        h_xx[i] = z2 / t2 - gpd_w_d2xi(x, z, t);
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
