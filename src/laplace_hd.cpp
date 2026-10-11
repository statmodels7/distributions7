#include <Rcpp.h>
#include "d7_par.h"
#include "ls_order_list.h"
using namespace Rcpp;

// Third/fourth-order derivatives of the Laplace log-density, almost everywhere:
// with r = y - mu, s = sign(r), a = |r|, the log-likelihood is
// l = -log(2b) - a/b, so every derivative is a monomial in a or s over a power
// of b. The kink at y = mu is the same one the gradient and Hessian carry;
// params_smooth already records it and the validators guard their references.

// [[Rcpp::export]]
List laplace_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector b,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_sigma(n), mu_sigma_sigma(n), sigma_sigma_sigma(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool b_is_scalar = (b.size() == 1);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double bb = b_is_scalar ? b[0] : b[i];
        double ib = 1.0 / bb, ib3 = ib * ib * ib;
        double r = y[i] - m;
        double s = (r > 0) - (r < 0);
        double u = std::abs(r) * ib;

        mu_mu_mu[i] = 0.0;
        mu_mu_sigma[i] = 0.0;
        mu_sigma_sigma[i] = 2.0 * s * ib3;
        sigma_sigma_sigma[i] = (6.0 * u - 2.0) * ib3;
    });

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_sigma") = mu_mu_sigma,
        Named("mu_sigma_sigma") = mu_sigma_sigma,
        Named("sigma_sigma_sigma") = sigma_sigma_sigma
    );
}

// [[Rcpp::export]]
List laplace_deriv3_expected_cpp(NumericVector y, NumericVector mu, NumericVector b,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_sigma(n), mu_sigma_sigma(n), sigma_sigma_sigma(n);
    bool b_is_scalar = (b.size() == 1);

    // E[s] = 0 and E[a] = b under the model, so only the pure-b component
    // survives: E[l_bbb] = -2/b^3 + 6/b^3 = 4/b^3.
    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double bb = b_is_scalar ? b[0] : b[i];
        double ib = 1.0 / bb;
        mu_mu_mu[i] = 0.0;
        mu_mu_sigma[i] = 0.0;
        mu_sigma_sigma[i] = 0.0;
        sigma_sigma_sigma[i] = 4.0 * ib * ib * ib;
    });

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_sigma") = mu_mu_sigma,
        Named("mu_sigma_sigma") = mu_sigma_sigma,
        Named("sigma_sigma_sigma") = sigma_sigma_sigma
    );
}

// [[Rcpp::export]]
List laplace_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector b,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_sigma(n), mu_mu_sigma_sigma(n), mu_sigma_sigma_sigma(n),
                  sigma_sigma_sigma_sigma(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool b_is_scalar = (b.size() == 1);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double bb = b_is_scalar ? b[0] : b[i];
        double ib = 1.0 / bb, ib2 = ib * ib, ib4 = ib2 * ib2;
        double r = y[i] - m;
        double s = (r > 0) - (r < 0);
        double u = std::abs(r) * ib;

        mu_mu_mu_mu[i] = 0.0;
        mu_mu_mu_sigma[i] = 0.0;
        mu_mu_sigma_sigma[i] = 0.0;
        mu_sigma_sigma_sigma[i] = -6.0 * s * ib4;
        sigma_sigma_sigma_sigma[i] = (6.0 - 24.0 * u) * ib4;
    });

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_sigma") = mu_mu_mu_sigma,
        Named("mu_mu_sigma_sigma") = mu_mu_sigma_sigma,
        Named("mu_sigma_sigma_sigma") = mu_sigma_sigma_sigma,
        Named("sigma_sigma_sigma_sigma") = sigma_sigma_sigma_sigma
    );
}

// [[Rcpp::export]]
List laplace_deriv4_expected_cpp(NumericVector y, NumericVector mu, NumericVector b,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_sigma(n), mu_mu_sigma_sigma(n), mu_sigma_sigma_sigma(n),
                  sigma_sigma_sigma_sigma(n);
    bool b_is_scalar = (b.size() == 1);

    // E[l_bbbb] = 6/b^4 - 24/b^4 = -18/b^4; every component carrying s is 0.
    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double bb = b_is_scalar ? b[0] : b[i];
        double ib = 1.0 / bb, ib2 = ib * ib;
        mu_mu_mu_mu[i] = 0.0;
        mu_mu_mu_sigma[i] = 0.0;
        mu_mu_sigma_sigma[i] = 0.0;
        mu_sigma_sigma_sigma[i] = 0.0;
        sigma_sigma_sigma_sigma[i] = -18.0 * ib2 * ib2;
    });

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_sigma") = mu_mu_mu_sigma,
        Named("mu_mu_sigma_sigma") = mu_mu_sigma_sigma,
        Named("mu_sigma_sigma_sigma") = mu_sigma_sigma_sigma,
        Named("sigma_sigma_sigma_sigma") = sigma_sigma_sigma_sigma
    );
}

// z^p g_i of the Laplace, g(z) = -log 2 - |z|: g1 = -sign(z), and every
// higher derivative is zero away from the kink
struct LaplaceTerm {
    double z;
    explicit LaplaceTerm(double zz) : z(zz) {}
    double operator()(int p, int i) const {
        if (i != 1) return 0.0;
        double sg = (z > 0) - (z < 0);
        return p == 0 ? -sg : -std::fabs(z);
    }
};

// [[Rcpp::export]]
List laplace_deriv5_cpp(NumericVector y, NumericVector mu, NumericVector b,
                        int threads = 1) {
    return d7::ls_order_list(5, y, mu, b, threads,
        [](double z) { return LaplaceTerm(z); });
}

// order 5 of l = log(lambda/2) - lambda|y - mu|: only the pure lambda
// component, 24/lambda^5, is nonzero away from the kink
// [[Rcpp::export]]
List laplace2_deriv5_cpp(NumericVector y, NumericVector mu, NumericVector lambda,
                         int threads = 1) {
    int n = y.size();
    NumericVector z0(n), z1(n), z2(n), z3(n), z4(n), o5(n);
    bool l_is_scalar = (lambda.size() == 1);
    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double il = 1.0 / (l_is_scalar ? lambda[0] : lambda[i]);
        double il2 = il * il;
        o5[i] = 24.0 * il2 * il2 * il;
    });
    return List::create(
        Named("mu_mu_mu_mu_mu") = z0,
        Named("mu_mu_mu_mu_lambda") = z1,
        Named("mu_mu_mu_lambda_lambda") = z2,
        Named("mu_mu_lambda_lambda_lambda") = z3,
        Named("mu_lambda_lambda_lambda_lambda") = z4,
        Named("lambda_lambda_lambda_lambda_lambda") = o5
    );
}
