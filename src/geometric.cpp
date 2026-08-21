#include <Rcpp.h>
#include "d7_par.h"
using namespace Rcpp;

// Geometric on {0, 1, 2, ...} in the MEAN parametrization. With success
// probability p = 1/(1+mu) the mass is p (1-p)^y, so
//   l = y log(mu) - (y + 1) log(1 + mu)
// and every order is the same shape twice, once at mu and once at 1 + mu:
//   l^(j) = (-1)^(j-1) (j-1)! [ y / mu^j - (y + 1) / (1 + mu)^j ].
// At j = 1 this collapses to (y - mu) / (mu (1 + mu)), the score of a
// one-parameter family in its mean over its variance. Under E[y] = mu,
//   E[l^(j)] = (-1)^(j-1) (j-1)! [ mu^(1-j) - (1 + mu)^(1-j) ],
// zero at j = 1 and -1/(mu(1+mu)) at j = 2.

static inline void geom_terms(double y, double m, double &a, double &b) {
    a = y;
    b = y + 1.0;
}

// [[Rcpp::export]]
List geometric_gradient_cpp(NumericVector y, NumericVector mu,
                        int threads = 1) {
    int n = y.size();
    NumericVector g(n);
    bool mu_is_scalar = (mu.size() == 1);

    d7::par_for(n, threads, d7::kMinTiny, [&](std::size_t i) {
        // LOCAL to the region: a scalar hoisted out of the loop and
        // written inside it is shared once the iterations are split
        const double m = mu_is_scalar ? mu[0] : mu[i];
        g[i] = (y[i] - m) / (m * (1.0 + m));
    });
    return List::create(Named("mu") = g);
}

// [[Rcpp::export]]
List geometric_hessian_cpp(NumericVector y, NumericVector mu,
                        int threads = 1) {
    int n = y.size();
    NumericVector h(n);
    bool mu_is_scalar = (mu.size() == 1);

    d7::par_for(n, threads, d7::kMinTiny, [&](std::size_t i) {
        // LOCAL to the region: a scalar hoisted out of the loop and
        // written inside it is shared once the iterations are split
        const double m = mu_is_scalar ? mu[0] : mu[i];
        double a, b;
        geom_terms(y[i], m, a, b);
        double om = 1.0 + m;
        h[i] = -(a / (m * m) - b / (om * om));
    });
    return List::create(Named("mu_mu") = h);
}

// [[Rcpp::export]]
List geometric_expected_hessian_cpp(NumericVector y, NumericVector mu,
                        int threads = 1) {
    int n = y.size();
    NumericVector h(n);
    bool mu_is_scalar = (mu.size() == 1);

    d7::par_for(n, threads, d7::kMinTiny, [&](std::size_t i) {
        // LOCAL to the region: a scalar hoisted out of the loop and
        // written inside it is shared once the iterations are split
        const double m = mu_is_scalar ? mu[0] : mu[i];
        h[i] = -1.0 / (m * (1.0 + m));
    });
    return List::create(Named("mu_mu") = h);
}

// [[Rcpp::export]]
List geometric_deriv3_cpp(NumericVector y, NumericVector mu,
                        int threads = 1) {
    int n = y.size();
    NumericVector d(n);
    bool mu_is_scalar = (mu.size() == 1);

    d7::par_for(n, threads, d7::kMinTiny, [&](std::size_t i) {
        // LOCAL to the region: a scalar hoisted out of the loop and
        // written inside it is shared once the iterations are split
        const double m = mu_is_scalar ? mu[0] : mu[i];
        double a, b;
        geom_terms(y[i], m, a, b);
        double om = 1.0 + m;
        d[i] = 2.0 * (a / (m * m * m) - b / (om * om * om));
    });
    return List::create(Named("mu_mu_mu") = d);
}

// [[Rcpp::export]]
List geometric_deriv3_expected_cpp(NumericVector y, NumericVector mu,
                        int threads = 1) {
    int n = y.size();
    NumericVector d(n);
    bool mu_is_scalar = (mu.size() == 1);

    d7::par_for(n, threads, d7::kMinTiny, [&](std::size_t i) {
        // LOCAL to the region: a scalar hoisted out of the loop and
        // written inside it is shared once the iterations are split
        const double m = mu_is_scalar ? mu[0] : mu[i];
        double om = 1.0 + m;
        d[i] = 2.0 * (1.0 / (m * m) - 1.0 / (om * om));
    });
    return List::create(Named("mu_mu_mu") = d);
}

// [[Rcpp::export]]
List geometric_deriv4_cpp(NumericVector y, NumericVector mu,
                        int threads = 1) {
    int n = y.size();
    NumericVector d(n);
    bool mu_is_scalar = (mu.size() == 1);

    d7::par_for(n, threads, d7::kMinTiny, [&](std::size_t i) {
        // LOCAL to the region: a scalar hoisted out of the loop and
        // written inside it is shared once the iterations are split
        const double m = mu_is_scalar ? mu[0] : mu[i];
        double a, b;
        geom_terms(y[i], m, a, b);
        double m2 = m * m, om = 1.0 + m, om2 = om * om;
        d[i] = -6.0 * (a / (m2 * m2) - b / (om2 * om2));
    });
    return List::create(Named("mu_mu_mu_mu") = d);
}

// [[Rcpp::export]]
List geometric_deriv4_expected_cpp(NumericVector y, NumericVector mu,
                        int threads = 1) {
    int n = y.size();
    NumericVector d(n);
    bool mu_is_scalar = (mu.size() == 1);

    d7::par_for(n, threads, d7::kMinTiny, [&](std::size_t i) {
        // LOCAL to the region: a scalar hoisted out of the loop and
        // written inside it is shared once the iterations are split
        const double m = mu_is_scalar ? mu[0] : mu[i];
        double om = 1.0 + m;
        d[i] = -6.0 * (1.0 / (m * m * m) - 1.0 / (om * om * om));
    });
    return List::create(Named("mu_mu_mu_mu") = d);
}
