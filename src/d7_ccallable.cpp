#include <Rcpp.h>
#include "psi_diff.h"
#include <R_ext/Rdynload.h>
#include <cstring>
#include <cmath>

// The scalar C entry points of the fast route piano_parallel.txt section 2a
// describes: the score and the second derivative of the log-density in ONE
// parameter at ONE observation, on the parameter scale, which is what a
// score-driven filter reads at every step of its recursion. A consumer
// resolves them once with R_GetCCallable and its loop then calls plain
// function pointers, touching no R API -- the precondition for any thread
// near that loop.
//
// The families covered are keyed by their S7 CLASS name, which is
// unambiguous where a distrib_name is shared across parametrizations, and
// an unknown name answers -1: the consumer keeps its R callbacks, so
// coverage is a speed property and never a correctness one. Every
// expression MIRRORS the family's own vector kernel line for line
// (gaussian.cpp, gamma1.cpp) -- the same operations in the same order, so
// the fast route is bit-identical to the callback route, which a twin test
// asserts. The remaining compiled families take the same few lines each
// when a measurement names them; the sixteen families whose derivatives
// live in vectorized R have no C body to point at and stay on the
// callbacks (piano_parallel.txt, section 3bis).

extern "C" {

int d7_scalar_id(const char* cls) {
    if (std::strcmp(cls, "Gaussian1Distrib") == 0) return 0;
    if (std::strcmp(cls, "Gamma1Distrib") == 0) return 1;
    return -1;
}

// k is the 0-based index of the parameter among the family's own, th the
// full parameter vector at this observation; out[0] the score component,
// out[1] the (k, k) second derivative, both on the parameter scale
void d7_score_curv(int id, int k, double y, const double* th, double* out) {
    switch (id) {
    case 0: {                       // gaussian1 (mu, sigma), as gaussian.cpp
        double m = th[0], sdev = th[1];
        double inv = 1.0 / sdev;
        double z = (y - m) / sdev;
        if (k == 0) {
            double inv2 = inv * inv;
            out[0] = z * inv;
            out[1] = -inv2;
        } else {
            double inv2 = inv * inv;
            out[0] = (z * z - 1.0) * inv;
            out[1] = (1.0 - 3.0 * z * z) * inv2;
        }
        break;
    }
    case 1: {                       // gamma1 (mu, phi), as gamma1.cpp
        double m = th[0], p = th[1];
        double s = 1.0 / p;
        double z = y / m;
        if (k == 0) {
            double m2 = m * m;
            out[0] = s * (z - 1.0) / m;
            out[1] = s * (1.0 - 2.0 * z) / m2;
        } else {
            // the same expressions gamma1_parts writes for f1 and f2; the
            // higher polygammas that struct also carries are not needed
            // here, and skipping them changes no computed value
            // see psi_diff.h; the R method's expression, written out
            double f1 = d7::psi_log_rest(s) + d7::psi_Ew2(z, z - 1.0);
            double f2 = 1.0 / s - R::trigamma(s);
            double s1 = -s * s, s2 = 2.0 * s * s * s;
            out[0] = f1 * (-s * s);
            out[1] = f2 * s1 * s1 + f1 * s2;
        }
        break;
    }
    default:
        out[0] = R_NaN; out[1] = R_NaN;
    }
}

// The expected information a filter driven by a scaled score reads at every
// step: out[0] the (k, k) expected second derivative E[l_kk], out[1] its
// derivative in the same parameter, both on the parameter scale. Each line
// mirrors the family's expected_hessian and dexpected kernels (gaussian.cpp,
// gamma1.cpp), so the scaled fast route is the scaled callback route.
void d7_info_dinfo(int id, int k, double y, const double* th, double* out) {
    switch (id) {
    case 0: {                       // gaussian1 (mu, sigma), as gaussian.cpp
        double sdev = th[1];
        double inv = 1.0 / sdev;
        double inv2 = inv * inv;
        if (k == 0) {
            out[0] = -inv2;
            out[1] = 0.0;
        } else {
            double inv3 = inv * inv * inv;
            out[0] = -2.0 * inv2;
            out[1] = 4.0 * inv3;
        }
        break;
    }
    case 1: {                       // gamma1 (mu, phi), as gamma1.cpp
        double m = th[0], p = th[1];
        double s = 1.0 / p;
        if (k == 0) {
            double m2 = m * m;
            double im = 1.0 / m, im2 = im * im;
            out[0] = -s / m2;
            out[1] = 2.0 * s * im2 * im;
        } else {
            double s1 = -s * s;
            double s2 = s * s, s3 = s2 * s, s4 = s2 * s2;
            double f2 = d7::psi1_rest(s), f3 = d7::psi2_rest(s);
            double q1 = f3 * s4 + 4.0 * f2 * s3;
            out[0] = f2 * s1 * s1;
            out[1] = -s2 * q1;
        }
        break;
    }
    default:
        out[0] = R_NaN; out[1] = R_NaN;
    }
}

} // extern "C"

// exposed to this package's own tests: the twin comparison against
// distrib_gradient()/distrib_hessian() lives where the formulas do
// [[Rcpp::export]]
Rcpp::List d7_scalar_probe(std::string cls, int k, Rcpp::NumericVector y,
                           Rcpp::NumericMatrix theta) {
    int id = d7_scalar_id(cls.c_str());
    int n = y.size(), np = theta.ncol();
    Rcpp::NumericVector g(n), h(n);
    std::vector<double> th(np);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < np; ++j) th[j] = theta(i, j);
        double out[2];
        d7_score_curv(id, k - 1, y[i], th.data(), out);
        g[i] = out[0]; h[i] = out[1];
    }
    return Rcpp::List::create(Rcpp::_["id"] = id, Rcpp::_["score"] = g,
                              Rcpp::_["curvature"] = h);
}

// the twin of d7_scalar_probe() for the expected information
// [[Rcpp::export]]
Rcpp::List d7_info_probe(std::string cls, int k, Rcpp::NumericVector y,
                         Rcpp::NumericMatrix theta) {
    int id = d7_scalar_id(cls.c_str());
    int n = y.size(), np = theta.ncol();
    Rcpp::NumericVector e(n), de(n);
    std::vector<double> th(np);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < np; ++j) th[j] = theta(i, j);
        double out[2];
        d7_info_dinfo(id, k - 1, y[i], th.data(), out);
        e[i] = out[0]; de[i] = out[1];
    }
    return Rcpp::List::create(Rcpp::_["id"] = id, Rcpp::_["expected"] = e,
                              Rcpp::_["dexpected"] = de);
}

// [[Rcpp::init]]
void d7_register_ccallable(DllInfo* dll) {
    R_RegisterCCallable("distributions7", "d7_scalar_id",
                        (DL_FUNC) d7_scalar_id);
    R_RegisterCCallable("distributions7", "d7_score_curv",
                        (DL_FUNC) d7_score_curv);
    R_RegisterCCallable("distributions7", "d7_info_dinfo",
                        (DL_FUNC) d7_info_dinfo);
}
