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

// Bernoulli, E = -1/q with q = m(1-m).
// [[Rcpp::export]]
List bernoulli_dexpected_cpp(NumericVector y, NumericVector mu, int order,
                             int threads = 1) {
    Par M(mu);
    return dexp_run(y.size(), threads, d7::kMinCheap, dexp_keys1("mu", order),
                    [&](std::size_t i, double* o) {
        double m = M[i], q = m * (1.0 - m), s = 1.0 - 2.0 * m, q2 = q * q;
        o[0] = order == 1 ? s / q2 : -2.0 / q2 - 2.0 * s * s / (q2 * q);
    });
}

// Binomial, the Bernoulli's times the size.
// [[Rcpp::export]]
List binomial_dexpected_cpp(NumericVector y, NumericVector mu,
                            NumericVector size, int order, int threads = 1) {
    Par M(mu), N(size);
    return dexp_run(y.size(), threads, d7::kMinCheap, dexp_keys1("mu", order),
                    [&](std::size_t i, double* o) {
        double m = M[i], q = m * (1.0 - m), s = 1.0 - 2.0 * m, q2 = q * q;
        o[0] = N[i] * (order == 1 ? s / q2 : -2.0 / q2 - 2.0 * s * s / (q2 * q));
    });
}

// Exponential, E = -1/m^2.
// [[Rcpp::export]]
List exponential_dexpected_cpp(NumericVector y, NumericVector mu, int order,
                               int threads = 1) {
    Par M(mu);
    return dexp_run(y.size(), threads, d7::kMinCheap, dexp_keys1("mu", order),
                    [&](std::size_t i, double* o) {
        double u = 1.0 / M[i], u3 = u * u * u;
        o[0] = order == 1 ? 2.0 * u3 : -6.0 * u3 * u;
    });
}

// Geometric, E = -1/q with q = m(1+m).
// [[Rcpp::export]]
List geometric_dexpected_cpp(NumericVector y, NumericVector mu, int order,
                             int threads = 1) {
    Par M(mu);
    return dexp_run(y.size(), threads, d7::kMinCheap, dexp_keys1("mu", order),
                    [&](std::size_t i, double* o) {
        double m = M[i], q = m * (1.0 + m), s = 1.0 + 2.0 * m, q2 = q * q;
        o[0] = order == 1 ? s / q2 : 2.0 / q2 - 2.0 * s * s / (q2 * q);
    });
}

// Chi-squared, E = -psi'(m/2)/4.
// [[Rcpp::export]]
List chisq_dexpected_cpp(NumericVector y, NumericVector mu, int order,
                         int threads = 1) {
    Par M(mu);
    return dexp_run(y.size(), threads, d7::kMinCostly, dexp_keys1("mu", order),
                    [&](std::size_t i, double* o) {
        double h = 0.5 * M[i];
        o[0] = order == 1 ? -R::psigamma(h, 2) / 8.0 : -R::psigamma(h, 3) / 16.0;
    });
}

// ---- location and scale, E_ab = k_ab / sigma^2 ----------------------------

namespace {
// d/ds (k/s^2) = -2k/s^3, d2/ds2 = 6k/s^4; nothing moves with the location.
template <typename Kfun>
List locscale_dexpected(NumericVector y, NumericVector sigma, int order,
                        int threads, const Kfun& kab) {
    Par S(sigma);
    double k[3];
    kab(k);   // E_mm, E_ss, E_ms, each as k / sigma^2
    return dexp_run(y.size(), threads, d7::kMinCheap,
                    dexp_keys2("mu", "sigma", order),
                    [&](std::size_t i, double* o) {
        double u = 1.0 / S[i], u2 = u * u, u3 = u2 * u, u4 = u2 * u2;
        for (int r = 0; r < 3; ++r) {
            if (order == 1) {
                o[2 * r] = 0.0;
                o[2 * r + 1] = -2.0 * k[r] * u3;
            } else {
                o[3 * r] = 0.0;
                o[3 * r + 1] = 6.0 * k[r] * u4;
                o[3 * r + 2] = 0.0;
            }
        }
    });
}
}  // namespace

// [[Rcpp::export]]
List cauchy_dexpected_cpp(NumericVector y, NumericVector mu,
                          NumericVector sigma, int order, int threads = 1) {
    return locscale_dexpected(y, sigma, order, threads, [](double* k) {
        k[0] = -0.5; k[1] = -0.5; k[2] = 0.0;
    });
}

// [[Rcpp::export]]
List logistic_dexpected_cpp(NumericVector y, NumericVector mu,
                            NumericVector sigma, int order, int threads = 1) {
    return locscale_dexpected(y, sigma, order, threads, [](double* k) {
        k[0] = -1.0 / 3.0; k[1] = -(3.0 + M_PI * M_PI) / 9.0; k[2] = 0.0;
    });
}

// [[Rcpp::export]]
List gumbel_dexpected_cpp(NumericVector y, NumericVector mu,
                          NumericVector sigma, int order, int threads = 1) {
    return locscale_dexpected(y, sigma, order, threads, [](double* k) {
        k[0] = -1.0; k[1] = -kGumbelC; k[2] = 1.0 - kEG;
    });
}

// ---- two parameters ------------------------------------------------------

// Gaussian by its variance (also the lognormal in (mu, sigma2)):
// E_mm = -1/v, E_vv = -1/(2 v^2), E_mv = 0.
// [[Rcpp::export]]
List gaussian2_dexpected_cpp(NumericVector y, NumericVector mu,
                             NumericVector sigma2, int order,
                             int threads = 1) {
    Par V(sigma2);
    return dexp_run(y.size(), threads, d7::kMinCheap,
                    dexp_keys2("mu", "sigma2", order),
                    [&](std::size_t i, double* o) {
        double u = 1.0 / V[i], u2 = u * u, u3 = u2 * u, u4 = u2 * u2;
        if (order == 1) {
            double w[6] = {0.0, u2, 0.0, u3, 0.0, 0.0};
            for (int k = 0; k < 6; ++k) o[k] = w[k];
        } else {
            double w[9] = {0.0, -2.0 * u3, 0.0, 0.0, -3.0 * u4, 0.0, 0.0, 0.0, 0.0};
            for (int k = 0; k < 9; ++k) o[k] = w[k];
        }
    });
}

// Gaussian by its precision: E_mm = -t, E_tt = -1/(2 t^2), E_mt = 0.
// [[Rcpp::export]]
List gaussian3_dexpected_cpp(NumericVector y, NumericVector mu,
                             NumericVector tau, int order, int threads = 1) {
    Par T(tau);
    return dexp_run(y.size(), threads, d7::kMinCheap,
                    dexp_keys2("mu", "tau", order),
                    [&](std::size_t i, double* o) {
        double u = 1.0 / T[i], u2 = u * u, u3 = u2 * u, u4 = u2 * u2;
        if (order == 1) {
            double w[6] = {0.0, -1.0, 0.0, u3, 0.0, 0.0};
            for (int k = 0; k < 6; ++k) o[k] = w[k];
        } else {
            double w[9] = {0.0, 0.0, 0.0, 0.0, -3.0 * u4, 0.0, 0.0, 0.0, 0.0};
            for (int k = 0; k < 9; ++k) o[k] = w[k];
        }
    });
}

// Inverse gaussian by its dispersion:
// E_mm = -1/(phi m^3), E_pp = -1/(2 phi^2), E_mp = 0.
// [[Rcpp::export]]
List invgauss1_dexpected_cpp(NumericVector y, NumericVector mu,
                             NumericVector phi, int order, int threads = 1) {
    Par M(mu), F(phi);
    return dexp_run(y.size(), threads, d7::kMinCheap,
                    dexp_keys2("mu", "phi", order),
                    [&](std::size_t i, double* o) {
        double m = M[i], p = F[i];
        double im = 1.0 / m, im3 = im * im * im, im4 = im3 * im, im5 = im4 * im;
        double ip = 1.0 / p, ip2 = ip * ip, ip3 = ip2 * ip, ip4 = ip2 * ip2;
        if (order == 1) {
            double w[6] = {3.0 * ip * im4, ip2 * im3, 0.0, ip3, 0.0, 0.0};
            for (int k = 0; k < 6; ++k) o[k] = w[k];
        } else {
            double w[9] = {-12.0 * ip * im5, -2.0 * ip3 * im3, -3.0 * ip2 * im4,
                           0.0, -3.0 * ip4, 0.0, 0.0, 0.0, 0.0};
            for (int k = 0; k < 9; ++k) o[k] = w[k];
        }
    });
}

// Inverse gaussian by its shape: E_mm = -L/m^3, E_LL = -1/(2 L^2), E_mL = 0.
// [[Rcpp::export]]
List invgauss2_dexpected_cpp(NumericVector y, NumericVector mu,
                             NumericVector lambda, int order, int threads = 1) {
    Par M(mu), La(lambda);
    return dexp_run(y.size(), threads, d7::kMinCheap,
                    dexp_keys2("mu", "lambda", order),
                    [&](std::size_t i, double* o) {
        double m = M[i], L = La[i];
        double im = 1.0 / m, im3 = im * im * im, im4 = im3 * im, im5 = im4 * im;
        double iL = 1.0 / L, iL3 = iL * iL * iL, iL4 = iL3 * iL;
        if (order == 1) {
            double w[6] = {3.0 * L * im4, -im3, 0.0, iL3, 0.0, 0.0};
            for (int k = 0; k < 6; ++k) o[k] = w[k];
        } else {
            double w[9] = {-12.0 * L * im5, 0.0, 3.0 * im4,
                           0.0, -3.0 * iL4, 0.0, 0.0, 0.0, 0.0};
            for (int k = 0; k < 9; ++k) o[k] = w[k];
        }
    });
}

// Gamma by its mean and variance. With a = m^2/v and g = a psi'(a) - 1,
// E_mm = -(1 + 4g)/v, E_vv = -a g/v^2, E_mv = 2 m g/v^2. g and its first two
// derivatives are read through the rests of psi_diff.h, so nothing cancels at
// a large shape, where g ~ 1/(2a).
// [[Rcpp::export]]
List gamma2_dexpected_cpp(NumericVector y, NumericVector mu,
                          NumericVector sigma2, int order, int threads = 1) {
    Par M(mu), V(sigma2);
    return dexp_run(y.size(), threads, d7::kMinCostly,
                    dexp_keys2("mu", "sigma2", order),
                    [&](std::size_t i, double* o) {
        double m = M[i], v = V[i], a = m * m / v;
        double r1 = d7::psi1_rest(a), r2 = d7::psi2_rest(a);
        double g = -a * r1, g1 = -r1 - a * r2;
        double iv = 1.0 / v, iv2 = iv * iv, iv3 = iv2 * iv, iv4 = iv2 * iv2;
        if (order == 1) {
            o[0] = -8.0 * m * g1 * iv2;
            o[1] = (1.0 + 4.0 * g + 4.0 * a * g1) * iv2;
            o[2] = -2.0 * m * (g + a * g1) * iv3;
            o[3] = a * (3.0 * g + a * g1) * iv3;
            o[4] = (2.0 * g + 4.0 * a * g1) * iv2;
            o[5] = -2.0 * m * (2.0 * g + a * g1) * iv3;
            return;
        }
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
// E_ab = psi'(s), s = a + b.
// [[Rcpp::export]]
List beta2_dexpected_cpp(NumericVector y, NumericVector alpha,
                         NumericVector beta, int order, int threads = 1) {
    Par A(alpha), B(beta);
    return dexp_run(y.size(), threads, d7::kMinCostly,
                    dexp_keys2("alpha", "beta", order),
                    [&](std::size_t i, double* o) {
        double a = A[i], b = B[i], k = order + 1;
        double ps = R::psigamma(a + b, k), pa = R::psigamma(a, k),
               pb = R::psigamma(b, k);
        if (order == 1) {
            double w[6] = {ps - pa, ps, ps, ps - pb, ps, ps};
            for (int j = 0; j < 6; ++j) o[j] = w[j];
        } else {
            double w[9] = {ps - pa, ps, ps, ps, ps - pb, ps, ps, ps, ps};
            for (int j = 0; j < 9; ++j) o[j] = w[j];
        }
    });
}

// Weibull by its scale and shape: E_mm = -s^2/m^2, E_ss = -C/s^2,
// E_ms = (1 - gamma)/m, C = (1 - gamma)^2 + pi^2/6.
// [[Rcpp::export]]
List weibull1_dexpected_cpp(NumericVector y, NumericVector mu,
                            NumericVector sigma, int order, int threads = 1) {
    Par M(mu), S(sigma);
    const double k = 1.0 - kEG;
    return dexp_run(y.size(), threads, d7::kMinCheap,
                    dexp_keys2("mu", "sigma", order),
                    [&](std::size_t i, double* o) {
        double m = M[i], s = S[i];
        double im = 1.0 / m, im2 = im * im, im3 = im2 * im, im4 = im2 * im2;
        double is = 1.0 / s, is3 = is * is * is, is4 = is3 * is;
        if (order == 1) {
            double w[6] = {2.0 * s * s * im3, -2.0 * s * im2,
                           0.0, 2.0 * kGumbelC * is3,
                           -k * im2, 0.0};
            for (int j = 0; j < 6; ++j) o[j] = w[j];
        } else {
            double w[9] = {-6.0 * s * s * im4, -2.0 * im2, 4.0 * s * im3,
                           0.0, -6.0 * kGumbelC * is4, 0.0,
                           2.0 * k * im3, 0.0, 0.0};
            for (int j = 0; j < 9; ++j) o[j] = w[j];
        }
    });
}

// ---- three parameters ----------------------------------------------------

namespace {

// pairs in hess_names() order: A_A, B_B, C_C, A_B, A_C, B_C
std::vector<std::string> dexp_keys3(const std::string& A, const std::string& B,
                                    const std::string& C, int order) {
    std::string par[3] = {A, B, C};
    std::string pr[6] = {A + "_" + A, B + "_" + B, C + "_" + C,
                         A + "_" + B, A + "_" + C, B + "_" + C};
    std::vector<std::string> k;
    for (int i = 0; i < 6; ++i) {
        if (order == 1) {
            for (int j = 0; j < 3; ++j) k.push_back(pr[i] + "_" + par[j]);
        } else {
            for (int j = 0; j < 6; ++j) k.push_back(pr[i] + "_" + pr[j]);
        }
    }
    return k;
}

// The first and second derivatives in nu of E[l_nu_nu] of the Student t,
//   E = [psi'((v+1)/2) - psi'(v/2)]/4 + (v+5)/(2 v (v+1)(v+3)),
// whose terms cancel from order v^-2 down to v^-5 and v^-6. Below v = 30 the
// direct form, above it the series, whose coefficients are exact integers
// obtained from the duplication identity psi'((v+1)/2) - psi'(v/2) =
// 4 psi'(v) - 2 psi'(v/2) and the Bernoulli expansion of psi'. The two agree
// best between 30 and 50, at 2.6e-13 to 5.5e-13.
void t_Enunu_d(double v, double& e1, double& e2) {
    if (v >= 30.0) {
        const double u = 1.0 / v, u2 = u * u, u5 = u2 * u2 * u;
        e1 = u5 * (14.0 + u * (-65.0 + u * (237.0 + u * (-833.0 + u * (2908.0 +
             u * (-9909.0 + u * (32795.0 + u * (-107393.0 + u * (354282.0 +
             u * (-1164917.0 + u * (3720073.0 + u * (-11670705.0 +
             u * 38263736.0))))))))))));
        e2 = u5 * u * (-70.0 + u * (390.0 + u * (-1659.0 + u * (6664.0 +
             u * (-26172.0 + u * (99090.0 + u * (-360745.0 + u * (1288716.0 +
             u * (-4605666.0 + u * (16308838.0 + u * (-55801095.0 +
             u * (186731280.0 + u * -650483512.0))))))))))));
        return;
    }
    double q = v * v * v + 4.0 * v * v + 3.0 * v;
    double q1 = 3.0 * v * v + 8.0 * v + 3.0, q2 = 6.0 * v + 8.0, num = v + 5.0;
    double r1 = (q - num * q1) / (2.0 * q * q);
    double r2 = (-num * q2 * q - 2.0 * q1 * (q - num * q1)) / (2.0 * q * q * q);
    e1 = (R::psigamma(0.5 * (v + 1.0), 2) - R::psigamma(0.5 * v, 2)) / 8.0 + r1;
    e2 = (R::psigamma(0.5 * (v + 1.0), 3) - R::psigamma(0.5 * v, 3)) / 16.0 + r2;
}

}  // namespace

// Student t by location, scale and degrees of freedom. With w1 = 1 + 1/v and
// w3 = 1 + 3/v, E_mm = -(w1/w3)/s^2, E_ss = -(2/w3)/s^2,
// E_sn = 2 v^-2/(w1 w3)/s, E_nn as above, E_ms = E_mn = 0. Every factor of v is
// written through 1/v, so nothing is Inf/Inf at the nu a log link reaches.
// Each component is s^b H(v): d_s = b G/s, d_v = s^b H', ss = b(b-1) G/s^2,
// sv = b s^(b-1) H', vv = s^b H''.
// [[Rcpp::export]]
List student_t1_dexpected_cpp(NumericVector y, NumericVector mu,
                              NumericVector sigma, NumericVector nu, int order,
                              int threads = 1) {
    Par S(sigma), V(nu);
    return dexp_run(y.size(), threads, d7::kMinCostly,
                    dexp_keys3("mu", "sigma", "nu", order),
                    [&](std::size_t i, double* o) {
        double s = S[i], iv = 1.0 / V[i], iv2 = iv * iv, iv3 = iv2 * iv;
        double w1 = 1.0 + iv, w3 = 1.0 + 3.0 * iv, Q = w1 * w3, t2 = 1.0 + 2.0 * iv;
        double H[6][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0},
                          {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        double b[6] = {-2.0, -2.0, 0.0, 0.0, 0.0, -1.0};
        H[0][0] = -w1 / w3; H[0][1] = -2.0 * iv2 / (w3 * w3);
        H[0][2] = 4.0 * iv3 / (w3 * w3 * w3);
        H[1][0] = -2.0 / w3; H[1][1] = -6.0 * iv2 / (w3 * w3);
        H[1][2] = 12.0 * iv3 / (w3 * w3 * w3);
        t_Enunu_d(V[i], H[2][1], H[2][2]);
        H[5][0] = 2.0 * iv2 / Q; H[5][1] = -4.0 * iv3 * t2 / (Q * Q);
        H[5][2] = -4.0 * iv2 * iv2 / (Q * Q) + 16.0 * iv2 * iv2 * t2 * t2 / (Q * Q * Q);
        for (int r = 0; r < 6; ++r) {
            double sb = std::pow(s, b[r]);
            double G = sb * H[r][0];
            if (order == 1) {
                o[3 * r] = 0.0;
                o[3 * r + 1] = b[r] * G / s;
                o[3 * r + 2] = sb * H[r][1];
            } else {
                double* w = o + 6 * r;
                w[0] = 0.0;                                  // mm
                w[1] = b[r] * (b[r] - 1.0) * G / (s * s);    // ss
                w[2] = sb * H[r][2];                         // nn
                w[3] = 0.0;                                  // ms
                w[4] = 0.0;                                  // mn
                w[5] = b[r] * sb / s * H[r][1];              // sn
            }
        }
    });
}

// Generalized gamma (Stacy) by scale a and shapes d, p, with k = d/p. Each
// expected component is a^al p^be F(k), so with k_d = 1/p and k_p = -k/p
//   d_a = al G/a, d_d = a^al p^(be-1) F', d_p = a^al p^(be-1) (be F - k F'),
//   aa = al(al-1) G/a^2, dd = a^al p^(be-2) F'',
//   pp = a^al p^(be-2) (be(be-1) F - 2(be-1) k F' + k^2 F''),
//   ad = al d_d / a, ap = al d_p / a, dp = a^al p^(be-2) ((be-1) F' - k F'').
// [[Rcpp::export]]
List gengamma1_dexpected_cpp(NumericVector y, NumericVector a, NumericVector d,
                             NumericVector p, int order, int threads = 1) {
    Par Aa(a), Dd(d), Pp(p);
    return dexp_run(y.size(), threads, d7::kMinCostly,
                    dexp_keys3("a", "d", "p", order),
                    [&](std::size_t i, double* o) {
        double av = Aa[i], dv = Dd[i], pv = Pp[i], k = dv / pv;
        double p0 = R::digamma(k), p1 = R::trigamma(k), p2 = R::psigamma(k, 2);
        double P0 = R::digamma(k + 1.0), P1 = R::trigamma(k + 1.0),
               P2 = R::psigamma(k + 1.0, 2);
        double p3 = 0.0, P3 = 0.0;
        if (order == 2) { p3 = R::psigamma(k, 3); P3 = R::psigamma(k + 1.0, 3); }
        // (al, be, F, F', F'') per pair: a_a, d_d, p_p, a_d, a_p, d_p
        double al[6] = {-2.0, 0.0, 0.0, -1.0, -1.0, 0.0};
        double be[6] = {2.0, -2.0, -2.0, 0.0, 0.0, -2.0};
        double F[6][3] = {
            {-k, -1.0, 0.0},
            {-p1, -p2, -p3},
            {-(1.0 + 2.0 * k * p0 + k * k * p1 + k * (P0 * P0 + P1)),
             -(2.0 * p0 + 4.0 * k * p1 + k * k * p2 + P0 * P0 + 2.0 * k * P0 * P1 +
               P1 + k * P2),
             -(6.0 * p1 + 6.0 * k * p2 + k * k * p3 + 4.0 * P0 * P1 +
               2.0 * k * P1 * P1 + 2.0 * k * P0 * P2 + 2.0 * P2 + k * P3)},
            {-1.0, 0.0, 0.0},
            {k * (1.0 + P0), 1.0 + P0 + k * P1, 2.0 * P1 + k * P2},
            {p0 + k * p1, 2.0 * p1 + k * p2, 3.0 * p2 + k * p3}};
        for (int r = 0; r < 6; ++r) {
            double aa = std::pow(av, al[r]);
            double G = aa * std::pow(pv, be[r]) * F[r][0];
            double pb1 = aa * std::pow(pv, be[r] - 1.0);
            double gd = pb1 * F[r][1];
            double gp = pb1 * (be[r] * F[r][0] - k * F[r][1]);
            if (order == 1) {
                o[3 * r] = al[r] * G / av;
                o[3 * r + 1] = gd;
                o[3 * r + 2] = gp;
            } else {
                double pb2 = aa * std::pow(pv, be[r] - 2.0);
                double* w = o + 6 * r;
                w[0] = al[r] * (al[r] - 1.0) * G / (av * av);
                w[1] = pb2 * F[r][2];
                w[2] = pb2 * (be[r] * (be[r] - 1.0) * F[r][0] -
                              2.0 * (be[r] - 1.0) * k * F[r][1] + k * k * F[r][2]);
                w[3] = al[r] * gd / av;
                w[4] = al[r] * gp / av;
                w[5] = pb2 * ((be[r] - 1.0) * F[r][1] - k * F[r][2]);
            }
        }
    });
}

// Generalized Pareto by scale and shape. With d = 1 + 2 xi, o = 1 + xi and
// F = 1/(d o), E_ss = -1/(d s^2), E_sx = -F/s, E_xx = -2F, and
// F' = -(3 + 4 xi) F^2, F'' = -4 F^2 + 2 (3 + 4 xi)^2 F^3. The information
// exists only for xi > -1/2, where it is NA, and so are its derivatives.
// [[Rcpp::export]]
List gpd_dexpected_cpp(NumericVector y, NumericVector sigma, NumericVector xi,
                       int order, int threads = 1) {
    Par S(sigma), X(xi);
    return dexp_run(y.size(), threads, d7::kMinCheap,
                    dexp_keys2("sigma", "xi", order),
                    [&](std::size_t i, double* o) {
        double s = S[i], x = X[i];
        int K = order == 1 ? 6 : 9;
        if (x <= -0.5) {
            for (int k = 0; k < K; ++k) o[k] = NA_REAL;
            return;
        }
        double d = 1.0 + 2.0 * x, F = 1.0 / (d * (1.0 + x)), c = 3.0 + 4.0 * x;
        double F1 = -c * F * F, F2 = -4.0 * F * F + 2.0 * c * c * F * F * F;
        double is = 1.0 / s, is2 = is * is, is3 = is2 * is, is4 = is2 * is2;
        double id = 1.0 / d, id2 = id * id, id3 = id2 * id;
        if (order == 1) {
            double w[6] = {2.0 * id * is3, 2.0 * id2 * is2,   // E_ss: d_s, d_x
                           0.0, -2.0 * F1,                     // E_xx
                           F * is2, -F1 * is};                 // E_sx
            for (int k = 0; k < 6; ++k) o[k] = w[k];
        } else {
            double w[9] = {-6.0 * id * is4, -8.0 * id3 * is2, -4.0 * id2 * is3,
                           0.0, -2.0 * F2, 0.0,
                           -2.0 * F * is3, -F2 * is, F1 * is2};
            for (int k = 0; k < 9; ++k) o[k] = w[k];
        }
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
// Returns E[l_ab] (three keys, as hess_names()), then the order-1 components,
// then at order 2 the order-2 ones, so one call serves both orders.
// [[Rcpp::export]]
List betabinom_shapes_dexpected_cpp(NumericVector y, NumericVector alpha,
                                    NumericVector beta, double size, int order,
                                    int threads = 1) {
    Par Av(alpha), Bv(beta);
    const int N = static_cast<int>(size);
    std::vector<std::string> keys = {"alpha_alpha", "beta_beta", "alpha_beta"};
    std::vector<std::string> k1 = dexp_keys2("alpha", "beta", 1);
    keys.insert(keys.end(), k1.begin(), k1.end());
    if (order == 2) {
        std::vector<std::string> k2 = dexp_keys2("alpha", "beta", 2);
        keys.insert(keys.end(), k2.begin(), k2.end());
    }
    const double fk[5] = {0.0, 1.0, -1.0, 2.0, -6.0};
    const int pu[3] = {0, 1, 0}, pv[3] = {0, 1, 1};
    return dexp_run(y.size(), threads, d7::kMinCostly, keys,
                    [&](std::size_t i, double* o) {
        const double a = Av[i], b = Bv[i];
        // cumulative power sums in a and in b, index m = number of terms
        std::vector<double> SA(4 * (N + 1), 0.0), SB(4 * (N + 1), 0.0);
        double SC[5] = {0.0, 0.0, 0.0, 0.0, 0.0}, lp = 0.0;
        for (int m = 0; m < N; ++m) {
            double ia = 1.0 / (a + m), ib = 1.0 / (b + m), ic = 1.0 / (a + b + m);
            double pa = 1.0, pb = 1.0, pc = 1.0;
            for (int k = 1; k <= 4; ++k) {
                pa *= ia; pb *= ib; pc *= ic;
                SA[(k - 1) * (N + 1) + m + 1] = SA[(k - 1) * (N + 1) + m] + pa;
                SB[(k - 1) * (N + 1) + m + 1] = SB[(k - 1) * (N + 1) + m] + pb;
                SC[k] += pc;
            }
            lp -= std::log1p(a / (b + m));
        }
        double acc[18];
        for (int k = 0; k < 18; ++k) acc[k] = 0.0;
        for (int j = 0; j <= N; ++j) {
            // D[na][nb]: the derivative with na differentiations in a and nb in b
            double D[5][5];
            for (int na = 0; na <= 4; ++na) {
                for (int nb = 0; na + nb <= 4; ++nb) {
                    int k = na + nb;
                    if (k == 0) { D[na][nb] = 0.0; continue; }
                    double c = -fk[k] * SC[k];
                    if (nb == 0) c += fk[k] * SA[(k - 1) * (N + 1) + j];
                    else if (na == 0) c += fk[k] * SB[(k - 1) * (N + 1) + N - j];
                    D[na][nb] = c;
                }
            }
            auto L = [&](int n0, int n1) { return D[n0][n1]; };
            double p = std::exp(lp);
            int w = 0;
            for (int r = 0; r < 3; ++r) {
                int u = pu[r], v = pv[r];
                int n0 = (u == 0) + (v == 0), n1 = 2 - n0;
                acc[w++] += p * L(n0, n1);
            }
            for (int r = 0; r < 3; ++r) {
                int n0 = (pu[r] == 0) + (pv[r] == 0), n1 = 2 - n0;
                double lab = L(n0, n1);
                for (int c = 0; c < 2; ++c) {
                    int c0 = n0 + (c == 0), c1 = n1 + (c == 1);
                    double lc = L(c == 0, c == 1);
                    acc[w++] += p * (L(c0, c1) + lab * lc);
                }
            }
            if (order == 2) {
                for (int r = 0; r < 3; ++r) {
                    int n0 = (pu[r] == 0) + (pv[r] == 0), n1 = 2 - n0;
                    double lab = L(n0, n1);
                    for (int s = 0; s < 3; ++s) {
                        int c = pu[s], d = pv[s];
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
            }
            if (j < N) {
                lp += std::log((N - j) / (j + 1.0)) +
                      std::log((j + a) / (N - j - 1.0 + b));
            }
        }
        const int K = order == 2 ? 18 : 9;
        for (int k = 0; k < K; ++k) o[k] = acc[k];
    });
}
