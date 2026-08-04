#include <Rcpp.h>
using namespace Rcpp;

// Third/fourth-order derivatives of the Laplace log-density, almost everywhere:
// with r = y - mu, s = sign(r), a = |r|, the log-likelihood is
// l = -log(2b) - a/b, so every derivative is a monomial in a or s over a power
// of b. The kink at y = mu is the same one the gradient and Hessian carry;
// params_smooth already records it and the validators guard their references.

// [[Rcpp::export]]
List laplace_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector b) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_b(n), mu_b_b(n), b_b_b(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool b_is_scalar = (b.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double bb = b_is_scalar ? b[0] : b[i];
        double b3 = bb * bb * bb, b4 = b3 * bb;
        double r = y[i] - m;
        double s = (r > 0) - (r < 0);
        double a = std::abs(r);

        mu_mu_mu[i] = 0.0;
        mu_mu_b[i] = 0.0;
        mu_b_b[i] = 2.0 * s / b3;
        b_b_b[i] = -2.0 / b3 + 6.0 * a / b4;
    }

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_b") = mu_mu_b,
        Named("mu_b_b") = mu_b_b,
        Named("b_b_b") = b_b_b
    );
}

// [[Rcpp::export]]
List laplace_deriv3_expected_cpp(NumericVector y, NumericVector mu, NumericVector b) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_b(n), mu_b_b(n), b_b_b(n);
    bool b_is_scalar = (b.size() == 1);

    // E[s] = 0 and E[a] = b under the model, so only the pure-b component
    // survives: E[l_bbb] = -2/b^3 + 6/b^3 = 4/b^3.
    for (int i = 0; i < n; i++) {
        double bb = b_is_scalar ? b[0] : b[i];
        double b3 = bb * bb * bb;
        mu_mu_mu[i] = 0.0;
        mu_mu_b[i] = 0.0;
        mu_b_b[i] = 0.0;
        b_b_b[i] = 4.0 / b3;
    }

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_b") = mu_mu_b,
        Named("mu_b_b") = mu_b_b,
        Named("b_b_b") = b_b_b
    );
}

// [[Rcpp::export]]
List laplace_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector b) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_b(n), mu_mu_b_b(n), mu_b_b_b(n),
                  b_b_b_b(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool b_is_scalar = (b.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double bb = b_is_scalar ? b[0] : b[i];
        double b4 = bb * bb * bb * bb, b5 = b4 * bb;
        double r = y[i] - m;
        double s = (r > 0) - (r < 0);
        double a = std::abs(r);

        mu_mu_mu_mu[i] = 0.0;
        mu_mu_mu_b[i] = 0.0;
        mu_mu_b_b[i] = 0.0;
        mu_b_b_b[i] = -6.0 * s / b4;
        b_b_b_b[i] = 6.0 / b4 - 24.0 * a / b5;
    }

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_b") = mu_mu_mu_b,
        Named("mu_mu_b_b") = mu_mu_b_b,
        Named("mu_b_b_b") = mu_b_b_b,
        Named("b_b_b_b") = b_b_b_b
    );
}

// [[Rcpp::export]]
List laplace_deriv4_expected_cpp(NumericVector y, NumericVector mu, NumericVector b) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_b(n), mu_mu_b_b(n), mu_b_b_b(n),
                  b_b_b_b(n);
    bool b_is_scalar = (b.size() == 1);

    // E[l_bbbb] = 6/b^4 - 24/b^4 = -18/b^4; every component carrying s is 0.
    for (int i = 0; i < n; i++) {
        double bb = b_is_scalar ? b[0] : b[i];
        double b4 = bb * bb * bb * bb;
        mu_mu_mu_mu[i] = 0.0;
        mu_mu_mu_b[i] = 0.0;
        mu_mu_b_b[i] = 0.0;
        mu_b_b_b[i] = 0.0;
        b_b_b_b[i] = -18.0 / b4;
    }

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_b") = mu_mu_mu_b,
        Named("mu_mu_b_b") = mu_mu_b_b,
        Named("mu_b_b_b") = mu_b_b_b,
        Named("b_b_b_b") = b_b_b_b
    );
}
