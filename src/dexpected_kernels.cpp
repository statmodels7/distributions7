// The first and second derivatives of the expected information in the
// parameters, for the families whose expected information is an elementary
// function of them. Each kernel returns, on the parameter scale, the list
// dexpected_names() (order 1) or d2expected_names() (order 2) enumerates; the
// link scale is carried across once, in R, by dexpected_link().
//
// Every component is an ordinary derivative of the family's written-out
// E[l_ab]. They were derived by hand and checked, before being written here,
// against one Richardson difference of the analytic order below and against
// the moment identities d_c E_ab = E[l_abc] + E[l_ab l_c] and its second-order
// counterpart through expectation().
//
// One body per kernel, decomposed over the observations by d7::par_for(), so
// the result is identical at any thread count.

#include <Rcpp.h>
#include <string>
#include <vector>
#include "d7_par.h"
#include "psi_diff.h"
using namespace Rcpp;

namespace {

// The keys, in the order the bodies below write their values:
// one parameter  : A_A_A  /  A_A_A_A
// two parameters : pairs (A_A, B_B, A_B) crossed with (A, B) at order 1 and
//                  with the pairs again at order 2, pair-major.
std::vector<std::string> dexp_keys1(const std::string& A, int order) {
    std::string p = A + "_" + A;
    return std::vector<std::string>(1, order == 1 ? p + "_" + A : p + "_" + p);
}

std::vector<std::string> dexp_keys2(const std::string& A, const std::string& B,
                                    int order) {
    std::string pr[3] = {A + "_" + A, B + "_" + B, A + "_" + B};
    std::vector<std::string> k;
    for (int i = 0; i < 3; ++i) {
        if (order == 1) {
            k.push_back(pr[i] + "_" + A);
            k.push_back(pr[i] + "_" + B);
        } else {
            for (int j = 0; j < 3; ++j) k.push_back(pr[i] + "_" + pr[j]);
        }
    }
    return k;
}

// Runs body(i, out) over the observations; body writes keys.size() doubles.
template <typename Body>
List dexp_run(int n, int threads, int threshold,
              const std::vector<std::string>& keys, const Body& body) {
    const int K = keys.size();
    std::vector<NumericVector> v(K);
    std::vector<double*> p(K);
    for (int k = 0; k < K; ++k) {
        v[k] = NumericVector(n);
        p[k] = v[k].begin();
    }
    d7::par_for(n, threads, threshold, [&](std::size_t i) {
        double o[36];
        body(i, o);
        for (int k = 0; k < K; ++k) p[k][i] = o[k];
    });
    List out(K);
    CharacterVector nm(K);
    for (int k = 0; k < K; ++k) {
        out[k] = v[k];
        nm[k] = keys[k];
    }
    out.attr("names") = nm;
    return out;
}

// A parameter read at observation i, recycled when it is a scalar.
struct Par {
    const double* x;
    bool s;
    explicit Par(const NumericVector& v) : x(v.begin()), s(v.size() == 1) {}
    double operator[](std::size_t i) const { return s ? x[0] : x[i]; }
};

const double kEG = 0.57721566490153286061;          // Euler-Mascheroni
const double kGumbelC = (1.0 - kEG) * (1.0 - kEG) + M_PI * M_PI / 6.0;

}  // namespace

// ---- one parameter -------------------------------------------------------
// One kernel per order throughout: dexpected1 returns the first derivatives
// of the expected information, dexpected2 the second.

// Bernoulli, E = -1/q with q = m(1-m).
// [[Rcpp::export]]
List bernoulli_dexpected1_cpp(NumericVector y, NumericVector mu, int threads = 1) {
    Par M(mu);
    return dexp_run(y.size(), threads, d7::kMinCheap, dexp_keys1("mu", 1),
                    [&](std::size_t i, double* o) {
        double m = M[i], q = m * (1.0 - m), s = 1.0 - 2.0 * m, q2 = q * q;
        o[0] = s / q2;
    });
}

// [[Rcpp::export]]
List bernoulli_dexpected2_cpp(NumericVector y, NumericVector mu, int threads = 1) {
    Par M(mu);
    return dexp_run(y.size(), threads, d7::kMinCheap, dexp_keys1("mu", 2),
                    [&](std::size_t i, double* o) {
        double m = M[i], q = m * (1.0 - m), s = 1.0 - 2.0 * m, q2 = q * q;
        o[0] = -2.0 / q2 - 2.0 * s * s / (q2 * q);
    });
}

// Binomial, the Bernoulli's times the size.
// [[Rcpp::export]]
List binomial_dexpected1_cpp(NumericVector y, NumericVector mu,
                             NumericVector size, int threads = 1) {
    Par M(mu), N(size);
    return dexp_run(y.size(), threads, d7::kMinCheap, dexp_keys1("mu", 1),
                    [&](std::size_t i, double* o) {
        double m = M[i], q = m * (1.0 - m), s = 1.0 - 2.0 * m, q2 = q * q;
        o[0] = N[i] * (s / q2);
    });
}

// [[Rcpp::export]]
List binomial_dexpected2_cpp(NumericVector y, NumericVector mu,
                             NumericVector size, int threads = 1) {
    Par M(mu), N(size);
    return dexp_run(y.size(), threads, d7::kMinCheap, dexp_keys1("mu", 2),
                    [&](std::size_t i, double* o) {
        double m = M[i], q = m * (1.0 - m), s = 1.0 - 2.0 * m, q2 = q * q;
        o[0] = N[i] * (-2.0 / q2 - 2.0 * s * s / (q2 * q));
    });
}

// Exponential, E = -1/m^2.
// [[Rcpp::export]]
List exponential_dexpected1_cpp(NumericVector y, NumericVector mu,
                                int threads = 1) {
    Par M(mu);
    return dexp_run(y.size(), threads, d7::kMinCheap, dexp_keys1("mu", 1),
                    [&](std::size_t i, double* o) {
        double u = 1.0 / M[i], u3 = u * u * u;
        o[0] = 2.0 * u3;
    });
}

// [[Rcpp::export]]
List exponential_dexpected2_cpp(NumericVector y, NumericVector mu,
                                int threads = 1) {
    Par M(mu);
    return dexp_run(y.size(), threads, d7::kMinCheap, dexp_keys1("mu", 2),
                    [&](std::size_t i, double* o) {
        double u = 1.0 / M[i], u3 = u * u * u;
        o[0] = -6.0 * u3 * u;
    });
}

// Geometric, E = -1/q with q = m(1+m).
// [[Rcpp::export]]
List geometric_dexpected1_cpp(NumericVector y, NumericVector mu,
                              int threads = 1) {
    Par M(mu);
    return dexp_run(y.size(), threads, d7::kMinCheap, dexp_keys1("mu", 1),
                    [&](std::size_t i, double* o) {
        double m = M[i], q = m * (1.0 + m), s = 1.0 + 2.0 * m, q2 = q * q;
        o[0] = s / q2;
    });
}

// [[Rcpp::export]]
List geometric_dexpected2_cpp(NumericVector y, NumericVector mu,
                              int threads = 1) {
    Par M(mu);
    return dexp_run(y.size(), threads, d7::kMinCheap, dexp_keys1("mu", 2),
                    [&](std::size_t i, double* o) {
        double m = M[i], q = m * (1.0 + m), s = 1.0 + 2.0 * m, q2 = q * q;
        o[0] = 2.0 / q2 - 2.0 * s * s / (q2 * q);
    });
}

// Chi-squared, E = -psi'(m/2)/4.
// [[Rcpp::export]]
List chisq_dexpected1_cpp(NumericVector y, NumericVector mu, int threads = 1) {
    Par M(mu);
    return dexp_run(y.size(), threads, d7::kMinCostly, dexp_keys1("mu", 1),
                    [&](std::size_t i, double* o) {
        o[0] = -R::psigamma(0.5 * M[i], 2) / 8.0;
    });
}

// [[Rcpp::export]]
List chisq_dexpected2_cpp(NumericVector y, NumericVector mu, int threads = 1) {
    Par M(mu);
    return dexp_run(y.size(), threads, d7::kMinCostly, dexp_keys1("mu", 2),
                    [&](std::size_t i, double* o) {
        o[0] = -R::psigamma(0.5 * M[i], 3) / 16.0;
    });
}

// ---- location and scale, E_ab = k_ab / sigma^2 ----------------------------

namespace {
// d/ds (k/s^2) = -2k/s^3, d2/ds2 = 6k/s^4; nothing moves with the location.
template <typename Kfun>
List locscale_dexpected1(NumericVector y, NumericVector sigma, int threads,
                         const Kfun& kab) {
    Par S(sigma);
    double k[3];
    kab(k);   // E_mm, E_ss, E_ms, each as k / sigma^2
    return dexp_run(y.size(), threads, d7::kMinCheap,
                    dexp_keys2("mu", "sigma", 1),
                    [&](std::size_t i, double* o) {
        double u = 1.0 / S[i], u3 = u * u * u;
        for (int r = 0; r < 3; ++r) {
            o[2 * r] = 0.0;
            o[2 * r + 1] = -2.0 * k[r] * u3;
        }
    });
}

template <typename Kfun>
List locscale_dexpected2(NumericVector y, NumericVector sigma, int threads,
                         const Kfun& kab) {
    Par S(sigma);
    double k[3];
    kab(k);
    return dexp_run(y.size(), threads, d7::kMinCheap,
                    dexp_keys2("mu", "sigma", 2),
                    [&](std::size_t i, double* o) {
        double u = 1.0 / S[i], u2 = u * u, u4 = u2 * u2;
        for (int r = 0; r < 3; ++r) {
            o[3 * r] = 0.0;
            o[3 * r + 1] = 6.0 * k[r] * u4;
            o[3 * r + 2] = 0.0;
        }
    });
}

inline void cauchy_k(double* k) { k[0] = -0.5; k[1] = -0.5; k[2] = 0.0; }
inline void logistic_k(double* k) {
    k[0] = -1.0 / 3.0; k[1] = -(3.0 + M_PI * M_PI) / 9.0; k[2] = 0.0;
}
inline void gumbel_k(double* k) { k[0] = -1.0; k[1] = -kGumbelC; k[2] = 1.0 - kEG; }
}  // namespace

// [[Rcpp::export]]
List cauchy_dexpected1_cpp(NumericVector y, NumericVector mu,
                           NumericVector sigma, int threads = 1) {
    return locscale_dexpected1(y, sigma, threads, cauchy_k);
}

// [[Rcpp::export]]
List cauchy_dexpected2_cpp(NumericVector y, NumericVector mu,
                           NumericVector sigma, int threads = 1) {
    return locscale_dexpected2(y, sigma, threads, cauchy_k);
}

// [[Rcpp::export]]
List logistic_dexpected1_cpp(NumericVector y, NumericVector mu,
                             NumericVector sigma, int threads = 1) {
    return locscale_dexpected1(y, sigma, threads, logistic_k);
}

// [[Rcpp::export]]
List logistic_dexpected2_cpp(NumericVector y, NumericVector mu,
                             NumericVector sigma, int threads = 1) {
    return locscale_dexpected2(y, sigma, threads, logistic_k);
}

// [[Rcpp::export]]
List gumbel_dexpected1_cpp(NumericVector y, NumericVector mu,
                           NumericVector sigma, int threads = 1) {
    return locscale_dexpected1(y, sigma, threads, gumbel_k);
}

// [[Rcpp::export]]
List gumbel_dexpected2_cpp(NumericVector y, NumericVector mu,
                           NumericVector sigma, int threads = 1) {
    return locscale_dexpected2(y, sigma, threads, gumbel_k);
}

// ---- two parameters ------------------------------------------------------

// Gaussian by its variance (also the lognormal in (mu, sigma2)):
// E_mm = -1/v, E_vv = -1/(2 v^2), E_mv = 0.
// [[Rcpp::export]]
List gaussian2_dexpected1_cpp(NumericVector y, NumericVector mu,
                              NumericVector sigma2, int threads = 1) {
    Par V(sigma2);
    return dexp_run(y.size(), threads, d7::kMinCheap,
                    dexp_keys2("mu", "sigma2", 1),
                    [&](std::size_t i, double* o) {
        double u = 1.0 / V[i], u2 = u * u, u3 = u2 * u;
        double w[6] = {0.0, u2, 0.0, u3, 0.0, 0.0};
        for (int k = 0; k < 6; ++k) o[k] = w[k];
    });
}

// [[Rcpp::export]]
List gaussian2_dexpected2_cpp(NumericVector y, NumericVector mu,
                              NumericVector sigma2, int threads = 1) {
    Par V(sigma2);
    return dexp_run(y.size(), threads, d7::kMinCheap,
                    dexp_keys2("mu", "sigma2", 2),
                    [&](std::size_t i, double* o) {
        double u = 1.0 / V[i], u2 = u * u, u3 = u2 * u, u4 = u2 * u2;
        double w[9] = {0.0, -2.0 * u3, 0.0, 0.0, -3.0 * u4, 0.0, 0.0, 0.0, 0.0};
        for (int k = 0; k < 9; ++k) o[k] = w[k];
    });
}

// Gaussian by its precision: E_mm = -t, E_tt = -1/(2 t^2), E_mt = 0.
// [[Rcpp::export]]
List gaussian3_dexpected1_cpp(NumericVector y, NumericVector mu,
                              NumericVector tau, int threads = 1) {
    Par T(tau);
    return dexp_run(y.size(), threads, d7::kMinCheap,
                    dexp_keys2("mu", "tau", 1),
                    [&](std::size_t i, double* o) {
        double u = 1.0 / T[i], u3 = u * u * u;
        double w[6] = {0.0, -1.0, 0.0, u3, 0.0, 0.0};
        for (int k = 0; k < 6; ++k) o[k] = w[k];
    });
}

// [[Rcpp::export]]
List gaussian3_dexpected2_cpp(NumericVector y, NumericVector mu,
                              NumericVector tau, int threads = 1) {
    Par T(tau);
    return dexp_run(y.size(), threads, d7::kMinCheap,
                    dexp_keys2("mu", "tau", 2),
                    [&](std::size_t i, double* o) {
        double u = 1.0 / T[i], u2 = u * u, u4 = u2 * u2;
        double w[9] = {0.0, 0.0, 0.0, 0.0, -3.0 * u4, 0.0, 0.0, 0.0, 0.0};
        for (int k = 0; k < 9; ++k) o[k] = w[k];
    });
}

// Inverse gaussian by its dispersion:
// E_mm = -1/(phi m^3), E_pp = -1/(2 phi^2), E_mp = 0.
// [[Rcpp::export]]
List invgauss1_dexpected1_cpp(NumericVector y, NumericVector mu,
                              NumericVector phi, int threads = 1) {
    Par M(mu), F(phi);
    return dexp_run(y.size(), threads, d7::kMinCheap,
                    dexp_keys2("mu", "phi", 1),
                    [&](std::size_t i, double* o) {
        double m = M[i], p = F[i];
        double im = 1.0 / m, im3 = im * im * im, im4 = im3 * im;
        double ip = 1.0 / p, ip2 = ip * ip, ip3 = ip2 * ip;
        double w[6] = {3.0 * ip * im4, ip2 * im3, 0.0, ip3, 0.0, 0.0};
        for (int k = 0; k < 6; ++k) o[k] = w[k];
    });
}

// [[Rcpp::export]]
List invgauss1_dexpected2_cpp(NumericVector y, NumericVector mu,
                              NumericVector phi, int threads = 1) {
    Par M(mu), F(phi);
    return dexp_run(y.size(), threads, d7::kMinCheap,
                    dexp_keys2("mu", "phi", 2),
                    [&](std::size_t i, double* o) {
        double m = M[i], p = F[i];
        double im = 1.0 / m, im3 = im * im * im, im4 = im3 * im, im5 = im4 * im;
        double ip = 1.0 / p, ip2 = ip * ip, ip3 = ip2 * ip, ip4 = ip2 * ip2;
        double w[9] = {-12.0 * ip * im5, -2.0 * ip3 * im3, -3.0 * ip2 * im4,
                       0.0, -3.0 * ip4, 0.0, 0.0, 0.0, 0.0};
        for (int k = 0; k < 9; ++k) o[k] = w[k];
    });
}

// Inverse gaussian by its shape: E_mm = -L/m^3, E_LL = -1/(2 L^2), E_mL = 0.
// [[Rcpp::export]]
List invgauss2_dexpected1_cpp(NumericVector y, NumericVector mu,
                              NumericVector lambda, int threads = 1) {
    Par M(mu), La(lambda);
    return dexp_run(y.size(), threads, d7::kMinCheap,
                    dexp_keys2("mu", "lambda", 1),
                    [&](std::size_t i, double* o) {
        double m = M[i], L = La[i];
        double im = 1.0 / m, im3 = im * im * im, im4 = im3 * im;
        double iL = 1.0 / L, iL3 = iL * iL * iL;
        double w[6] = {3.0 * L * im4, -im3, 0.0, iL3, 0.0, 0.0};
        for (int k = 0; k < 6; ++k) o[k] = w[k];
    });
}

// [[Rcpp::export]]
List invgauss2_dexpected2_cpp(NumericVector y, NumericVector mu,
                              NumericVector lambda, int threads = 1) {
    Par M(mu), La(lambda);
    return dexp_run(y.size(), threads, d7::kMinCheap,
                    dexp_keys2("mu", "lambda", 2),
                    [&](std::size_t i, double* o) {
        double m = M[i], L = La[i];
        double im = 1.0 / m, im3 = im * im * im, im4 = im3 * im, im5 = im4 * im;
        double iL = 1.0 / L, iL3 = iL * iL * iL, iL4 = iL3 * iL;
        double w[9] = {-12.0 * L * im5, 0.0, 3.0 * im4,
                       0.0, -3.0 * iL4, 0.0, 0.0, 0.0, 0.0};
        for (int k = 0; k < 9; ++k) o[k] = w[k];
    });
}

// Gamma by its mean and variance. With a = m^2/v and g = a psi'(a) - 1,
// E_mm = -(1 + 4g)/v, E_vv = -a g/v^2, E_mv = 2 m g/v^2. g and its first two
// derivatives are read through the rests of psi_diff.h, so nothing cancels at
// a large shape, where g ~ 1/(2a).
// [[Rcpp::export]]
List gamma2_dexpected1_cpp(NumericVector y, NumericVector mu,
                           NumericVector sigma2, int threads = 1) {
    Par M(mu), V(sigma2);
    return dexp_run(y.size(), threads, d7::kMinCostly,
                    dexp_keys2("mu", "sigma2", 1),
                    [&](std::size_t i, double* o) {
        double m = M[i], v = V[i], a = m * m / v;
        double r1 = d7::psi1_rest(a), r2 = d7::psi2_rest(a);
        double g = -a * r1, g1 = -r1 - a * r2;
        double iv = 1.0 / v, iv2 = iv * iv, iv3 = iv2 * iv;
        o[0] = -8.0 * m * g1 * iv2;
        o[1] = (1.0 + 4.0 * g + 4.0 * a * g1) * iv2;
        o[2] = -2.0 * m * (g + a * g1) * iv3;
        o[3] = a * (3.0 * g + a * g1) * iv3;
        o[4] = (2.0 * g + 4.0 * a * g1) * iv2;
        o[5] = -2.0 * m * (2.0 * g + a * g1) * iv3;
    });
}

// [[Rcpp::export]]
List gamma2_dexpected2_cpp(NumericVector y, NumericVector mu,
                           NumericVector sigma2, int threads = 1) {
    Par M(mu), V(sigma2);
    return dexp_run(y.size(), threads, d7::kMinCostly,
                    dexp_keys2("mu", "sigma2", 2),
                    [&](std::size_t i, double* o) {
        double m = M[i], v = V[i], a = m * m / v;
        double r1 = d7::psi1_rest(a), r2 = d7::psi2_rest(a);
        double g = -a * r1, g1 = -r1 - a * r2;
        double iv = 1.0 / v, iv2 = iv * iv, iv3 = iv2 * iv, iv4 = iv2 * iv2;
        double r3 = d7::psi3_rest(a);
        double g2 = -2.0 * r2 - a * r3, a2 = a * a;
        o[0] = -(8.0 * g1 + 16.0 * a * g2) * iv2;
        o[1] = -(2.0 + 8.0 * g + 16.0 * a * g1 + 4.0 * a2 * g2) * iv3;
        o[2] = 8.0 * m * (a * g2 + 2.0 * g1) * iv3;
        o[3] = -(2.0 * g + 10.0 * a * g1 + 4.0 * a2 * g2) * iv3;
        o[4] = -a * (12.0 * g + 8.0 * a * g1 + a2 * g2) * iv4;
        o[5] = 2.0 * m * (3.0 * g + 5.0 * a * g1 + a2 * g2) * iv4;
        o[6] = 2.0 * m * (6.0 * g1 + 4.0 * a * g2) * iv3;
        o[7] = 2.0 * m * (6.0 * g + 6.0 * a * g1 + a2 * g2) * iv4;
        o[8] = -(4.0 * g + 14.0 * a * g1 + 4.0 * a2 * g2) * iv3;
    });
}

// Beta by its two shapes: E_aa = psi'(s) - psi'(a), E_bb = psi'(s) - psi'(b),
// E_ab = psi'(s), s = a + b; derivative k reads psi^(k+1).
// [[Rcpp::export]]
List beta2_dexpected1_cpp(NumericVector y, NumericVector alpha,
                          NumericVector beta, int threads = 1) {
    Par A(alpha), B(beta);
    return dexp_run(y.size(), threads, d7::kMinCostly,
                    dexp_keys2("alpha", "beta", 1),
                    [&](std::size_t i, double* o) {
        double a = A[i], b = B[i];
        double ps = R::psigamma(a + b, 2), pa = R::psigamma(a, 2),
               pb = R::psigamma(b, 2);
        double w[6] = {ps - pa, ps, ps, ps - pb, ps, ps};
        for (int j = 0; j < 6; ++j) o[j] = w[j];
    });
}

// [[Rcpp::export]]
List beta2_dexpected2_cpp(NumericVector y, NumericVector alpha,
                          NumericVector beta, int threads = 1) {
    Par A(alpha), B(beta);
    return dexp_run(y.size(), threads, d7::kMinCostly,
                    dexp_keys2("alpha", "beta", 2),
                    [&](std::size_t i, double* o) {
        double a = A[i], b = B[i];
        double ps = R::psigamma(a + b, 3), pa = R::psigamma(a, 3),
               pb = R::psigamma(b, 3);
        double w[9] = {ps - pa, ps, ps, ps, ps - pb, ps, ps, ps, ps};
        for (int j = 0; j < 9; ++j) o[j] = w[j];
    });
}

// Weibull by its scale and shape: E_mm = -s^2/m^2, E_ss = -C/s^2,
// E_ms = (1 - gamma)/m, C = (1 - gamma)^2 + pi^2/6.
// [[Rcpp::export]]
List weibull1_dexpected1_cpp(NumericVector y, NumericVector mu,
                             NumericVector sigma, int threads = 1) {
    Par M(mu), S(sigma);
    const double k = 1.0 - kEG;
    return dexp_run(y.size(), threads, d7::kMinCheap,
                    dexp_keys2("mu", "sigma", 1),
                    [&](std::size_t i, double* o) {
        double m = M[i], s = S[i];
        double im = 1.0 / m, im2 = im * im, im3 = im2 * im;
        double is = 1.0 / s, is3 = is * is * is;
        double w[6] = {2.0 * s * s * im3, -2.0 * s * im2,
                       0.0, 2.0 * kGumbelC * is3,
                       -k * im2, 0.0};
        for (int j = 0; j < 6; ++j) o[j] = w[j];
    });
}

// [[Rcpp::export]]
List weibull1_dexpected2_cpp(NumericVector y, NumericVector mu,
                             NumericVector sigma, int threads = 1) {
    Par M(mu), S(sigma);
    const double k = 1.0 - kEG;
    return dexp_run(y.size(), threads, d7::kMinCheap,
                    dexp_keys2("mu", "sigma", 2),
                    [&](std::size_t i, double* o) {
        double m = M[i], s = S[i];
        double im = 1.0 / m, im2 = im * im, im3 = im2 * im, im4 = im2 * im2;
        double is = 1.0 / s, is3 = is * is * is, is4 = is3 * is;
        double w[9] = {-6.0 * s * s * im4, -2.0 * im2, 4.0 * s * im3,
                       0.0, -6.0 * kGumbelC * is4, 0.0,
                       2.0 * k * im3, 0.0, 0.0};
        for (int j = 0; j < 9; ++j) o[j] = w[j];
    });
}

// Generalized Pareto by scale and shape. With d = 1 + 2 xi, o = 1 + xi and
// F = 1/(d o), E_ss = -1/(d s^2), E_sx = -F/s, E_xx = -2F, and
// F' = -(3 + 4 xi) F^2, F'' = -4 F^2 + 2 (3 + 4 xi)^2 F^3. The information
// exists only for xi > -1/2, where it is NA, and so are its derivatives.
// [[Rcpp::export]]
List gpd_dexpected1_cpp(NumericVector y, NumericVector sigma, NumericVector xi,
                        int threads = 1) {
    Par S(sigma), X(xi);
    return dexp_run(y.size(), threads, d7::kMinCheap,
                    dexp_keys2("sigma", "xi", 1),
                    [&](std::size_t i, double* o) {
        double s = S[i], x = X[i];
        if (x <= -0.5) {
            for (int k = 0; k < 6; ++k) o[k] = NA_REAL;
            return;
        }
        double d = 1.0 + 2.0 * x, F = 1.0 / (d * (1.0 + x)), c = 3.0 + 4.0 * x;
        double F1 = -c * F * F;
        double is = 1.0 / s, is2 = is * is, is3 = is2 * is;
        double id = 1.0 / d, id2 = id * id;
        double w[6] = {2.0 * id * is3, 2.0 * id2 * is2,   // E_ss: d_s, d_x
                       0.0, -2.0 * F1,                     // E_xx
                       F * is2, -F1 * is};                 // E_sx
        for (int k = 0; k < 6; ++k) o[k] = w[k];
    });
}

// [[Rcpp::export]]
List gpd_dexpected2_cpp(NumericVector y, NumericVector sigma, NumericVector xi,
                        int threads = 1) {
    Par S(sigma), X(xi);
    return dexp_run(y.size(), threads, d7::kMinCheap,
                    dexp_keys2("sigma", "xi", 2),
                    [&](std::size_t i, double* o) {
        double s = S[i], x = X[i];
        if (x <= -0.5) {
            for (int k = 0; k < 9; ++k) o[k] = NA_REAL;
            return;
        }
        double d = 1.0 + 2.0 * x, F = 1.0 / (d * (1.0 + x)), c = 3.0 + 4.0 * x;
        double F1 = -c * F * F, F2 = -4.0 * F * F + 2.0 * c * c * F * F * F;
        double is = 1.0 / s, is2 = is * is, is3 = is2 * is, is4 = is2 * is2;
        double id = 1.0 / d, id2 = id * id, id3 = id2 * id;
        double w[9] = {-6.0 * id * is4, -8.0 * id3 * is2, -4.0 * id2 * is3,
                       0.0, -2.0 * F2, 0.0,
                       -2.0 * F * is3, -F2 * is, F1 * is2};
        for (int k = 0; k < 9; ++k) o[k] = w[k];
    });
}

// ---- a finite support ----------------------------------------------------

// Beta-binomial in its shapes (a, b) with size N. Every expectation is a
// finite sum over 0..N against the mass, and every derivative of the log-mass
// is a sum of negative powers: with S_k(x, m) = sum_{i<m} (x+i)^-k and
// f_k = (-1)^(k-1) (k-1)!, the k-th derivative in a alone at y is
// f_k [S_k(a, y) - S_k(a+b, N)], in b alone f_k [S_k(b, N-y) - S_k(a+b, N)],
// and any mixed one -f_k S_k(a+b, N) -- the polygamma differences at integer
// shifts, accumulated exactly beside the mass, so no polygamma is called and
// none cancels. The mass starts from P(Y=0) = prod_{i<N} (b+i)/(a+b+i) on the
// log scale and follows the ratio of consecutive masses.
// One kernel per order. The first returns E[l_ab] (three keys, as
// hess_names()) and the first derivatives, which read the log-mass's
// derivatives to total order three; the second returns the second
// derivatives alone, which read them to order four.
namespace {

const double kBBfk[5] = {0.0, 1.0, -1.0, 2.0, -6.0};
const int kBBpu[3] = {0, 1, 0}, kBBpv[3] = {0, 1, 1};

// the cumulative power sums to order K (3 or 4) and log P(Y = 0)
struct BBSums {
    int N, K;
    std::vector<double> SA, SB;
    double SC[5];
    double lp0;
    BBSums(double a, double b, int N_, int K_)
        : N(N_), K(K_), SA(K_ * (N_ + 1), 0.0), SB(K_ * (N_ + 1), 0.0),
          lp0(0.0) {
        for (int k = 0; k < 5; ++k) SC[k] = 0.0;
        for (int m = 0; m < N; ++m) {
            double ia = 1.0 / (a + m), ib = 1.0 / (b + m), ic = 1.0 / (a + b + m);
            double pa = 1.0, pb = 1.0, pc = 1.0;
            for (int k = 1; k <= K; ++k) {
                pa *= ia; pb *= ib; pc *= ic;
                SA[(k - 1) * (N + 1) + m + 1] = SA[(k - 1) * (N + 1) + m] + pa;
                SB[(k - 1) * (N + 1) + m + 1] = SB[(k - 1) * (N + 1) + m] + pb;
                SC[k] += pc;
            }
            lp0 -= std::log1p(a / (b + m));
        }
    }
    // the derivatives at y = j with na differentiations in a and nb in b,
    // na + nb <= K, into D
    void at(int j, double D[5][5]) const {
        for (int na = 0; na <= K; ++na) {
            for (int nb = 0; na + nb <= K; ++nb) {
                int k = na + nb;
                if (k == 0) { D[na][nb] = 0.0; continue; }
                double c = -kBBfk[k] * SC[k];
                if (nb == 0) c += kBBfk[k] * SA[(k - 1) * (N + 1) + j];
                else if (na == 0) c += kBBfk[k] * SB[(k - 1) * (N + 1) + N - j];
                D[na][nb] = c;
            }
        }
    }
};

inline double bb_next_lp(double lp, int j, int N, double a, double b) {
    return lp + std::log((N - j) / (j + 1.0)) +
        std::log((j + a) / (N - j - 1.0 + b));
}

}  // namespace

// [[Rcpp::export]]
List betabinom_shapes_dexpected1_cpp(NumericVector y, NumericVector alpha,
                                     NumericVector beta, double size,
                                     int threads = 1) {
    Par Av(alpha), Bv(beta);
    const int N = static_cast<int>(size);
    std::vector<std::string> keys = {"alpha_alpha", "beta_beta", "alpha_beta"};
    std::vector<std::string> k1 = dexp_keys2("alpha", "beta", 1);
    keys.insert(keys.end(), k1.begin(), k1.end());
    return dexp_run(y.size(), threads, d7::kMinCostly, keys,
                    [&](std::size_t i, double* o) {
        const double a = Av[i], b = Bv[i];
        const BBSums S(a, b, N, 3);
        double lp = S.lp0;
        double acc[9];
        for (int k = 0; k < 9; ++k) acc[k] = 0.0;
        for (int j = 0; j <= N; ++j) {
            double D[5][5];
            S.at(j, D);
            double p = std::exp(lp);
            int w = 0;
            for (int r = 0; r < 3; ++r) {
                int n0 = (kBBpu[r] == 0) + (kBBpv[r] == 0), n1 = 2 - n0;
                acc[w++] += p * D[n0][n1];
            }
            for (int r = 0; r < 3; ++r) {
                int n0 = (kBBpu[r] == 0) + (kBBpv[r] == 0), n1 = 2 - n0;
                double lab = D[n0][n1];
                for (int c = 0; c < 2; ++c) {
                    int c0 = n0 + (c == 0), c1 = n1 + (c == 1);
                    double lc = D[c == 0][c == 1];
                    acc[w++] += p * (D[c0][c1] + lab * lc);
                }
            }
            if (j < N) lp = bb_next_lp(lp, j, N, a, b);
        }
        for (int k = 0; k < 9; ++k) o[k] = acc[k];
    });
}

// [[Rcpp::export]]
List betabinom_shapes_dexpected2_cpp(NumericVector y, NumericVector alpha,
                                     NumericVector beta, double size,
                                     int threads = 1) {
    Par Av(alpha), Bv(beta);
    const int N = static_cast<int>(size);
    return dexp_run(y.size(), threads, d7::kMinCostly,
                    dexp_keys2("alpha", "beta", 2),
                    [&](std::size_t i, double* o) {
        const double a = Av[i], b = Bv[i];
        const BBSums S(a, b, N, 4);
        double lp = S.lp0;
        double acc[9];
        for (int k = 0; k < 9; ++k) acc[k] = 0.0;
        for (int j = 0; j <= N; ++j) {
            double D[5][5];
            S.at(j, D);
            auto L = [&](int n0, int n1) { return D[n0][n1]; };
            double p = std::exp(lp);
            int w = 0;
            for (int r = 0; r < 3; ++r) {
                int n0 = (kBBpu[r] == 0) + (kBBpv[r] == 0), n1 = 2 - n0;
                double lab = L(n0, n1);
                for (int s = 0; s < 3; ++s) {
                    int c = kBBpu[s], d = kBBpv[s];
                    double lc = L(c == 0, c == 1), ld = L(d == 0, d == 1);
                    double lcd = L((c == 0) + (d == 0), (c == 1) + (d == 1));
                    double labc = L(n0 + (c == 0), n1 + (c == 1));
                    double labd = L(n0 + (d == 0), n1 + (d == 1));
                    double labcd = L(n0 + (c == 0) + (d == 0),
                                     n1 + (c == 1) + (d == 1));
                    acc[w++] += p * (labcd + labd * lc + labc * ld +
                                     lab * lcd + lab * lc * ld);
                }
            }
            if (j < N) lp = bb_next_lp(lp, j, N, a, b);
        }
        for (int k = 0; k < 9; ++k) o[k] = acc[k];
    });
}
