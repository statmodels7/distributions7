#include <Rcpp.h>
using namespace Rcpp;

// Third/fourth-order derivatives of the Gamma log-density (mean/variance
// parameterization, alpha = mu^2/sigma2), transcribed from the Wolfram output.
// La = digamma(alpha), Ta = trigamma(alpha), P2 = psigamma(alpha,2),
// P3 = psigamma(alpha,3), LL = log(mu/sigma2) + log(y), DL = La - LL.

// [[Rcpp::export]]
List gamma_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector sigma2) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_sigma2(n), mu_sigma2_sigma2(n), sigma2_sigma2_sigma2(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool s_is_scalar = (sigma2.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s2 = s_is_scalar ? sigma2[0] : sigma2[i];
        double m2 = m * m, m4 = m2 * m2, m5 = m4 * m;
        double s2_2 = s2 * s2, s2_3 = s2_2 * s2, s2_4 = s2_2 * s2_2, s2_5 = s2_4 * s2, s2_6 = s2_3 * s2_3;
        double a = m2 / s2;
        double La = R::digamma(a), Ta = R::trigamma(a), P2 = R::psigamma(a, 2.0);
        double DL = La - (std::log(m / s2) + std::log(y[i]));

        mu_mu_mu[i] = (-8.0 * m4 * P2 + 2.0 * s2 * (s2 - 6.0 * m2 * Ta)) / (m * s2_3);
        mu_mu_sigma2[i] = (-5.0 * s2_2 + 4.0 * m4 * P2 + 2.0 * s2 * (s2 * DL + 5.0 * m2 * Ta)) / s2_4;
        mu_sigma2_sigma2[i] = -2.0 * (s2_2 * (-4.0 * m + y[i]) + m5 * P2 + 2.0 * m * s2 * (s2 * DL + 2.0 * m2 * Ta)) / s2_5;
        sigma2_sigma2_sigma2[i] = (m * (s2_2 * (-11.0 * m + 6.0 * y[i]) + m5 * P2 + 6.0 * m * s2 * (s2 * DL + m2 * Ta))) / s2_6;
    }

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_sigma2") = mu_mu_sigma2,
        Named("mu_sigma2_sigma2") = mu_sigma2_sigma2,
        Named("sigma2_sigma2_sigma2") = sigma2_sigma2_sigma2
    );
}

// [[Rcpp::export]]
List gamma_deriv3_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma2) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_sigma2(n), mu_sigma2_sigma2(n), sigma2_sigma2_sigma2(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool s_is_scalar = (sigma2.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s2 = s_is_scalar ? sigma2[0] : sigma2[i];
        double m2 = m * m, m3 = m2 * m, m4 = m2 * m2, m5 = m4 * m;
        double s2_2 = s2 * s2, s2_3 = s2_2 * s2, s2_4 = s2_2 * s2_2, s2_5 = s2_4 * s2, s2_6 = s2_3 * s2_3;
        double a = m2 / s2;
        double Ta = R::trigamma(a), P2 = R::psigamma(a, 2.0);

        mu_mu_mu[i] = 2.0 * (-4.0 * m4 * P2 + s2 * (s2 - 6.0 * m2 * Ta)) / (m * s2_3);
        mu_mu_sigma2[i] = (4.0 * m4 * P2 - 5.0 * s2 * (s2 - 2.0 * m2 * Ta)) / s2_4;
        mu_sigma2_sigma2[i] = (6.0 * m * s2_2 - 2.0 * m5 * P2 - 8.0 * m3 * s2 * Ta) / s2_5;
        sigma2_sigma2_sigma2[i] = (m2 * (m4 * P2 + s2 * (-5.0 * s2 + 6.0 * m2 * Ta))) / s2_6;
    }

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_sigma2") = mu_mu_sigma2,
        Named("mu_sigma2_sigma2") = mu_sigma2_sigma2,
        Named("sigma2_sigma2_sigma2") = sigma2_sigma2_sigma2
    );
}

// [[Rcpp::export]]
List gamma_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector sigma2) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_sigma2(n), mu_mu_sigma2_sigma2(n),
                  mu_sigma2_sigma2_sigma2(n), sigma2_sigma2_sigma2_sigma2(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool s_is_scalar = (sigma2.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s2 = s_is_scalar ? sigma2[0] : sigma2[i];
        double m2 = m * m, m3 = m2 * m, m4 = m2 * m2, m5 = m4 * m, m6 = m4 * m2, m7 = m6 * m;
        double s2_2 = s2 * s2, s2_3 = s2_2 * s2, s2_4 = s2_2 * s2_2,
               s2_5 = s2_4 * s2, s2_6 = s2_3 * s2_3, s2_7 = s2_6 * s2, s2_8 = s2_4 * s2_4;
        double a = m2 / s2;
        double Ta = R::trigamma(a), P2 = R::psigamma(a, 2.0), P3 = R::psigamma(a, 3.0);
        double DL = R::digamma(a) - (std::log(m / s2) + std::log(y[i]));

        mu_mu_mu_mu[i] = -2.0 * (24.0 * m4 * s2 * P2 + 8.0 * m6 * P3 + s2_2 * (s2 + 6.0 * m2 * Ta)) / (m2 * s2_4);
        mu_mu_mu_sigma2[i] = (36.0 * m4 * s2 * P2 + 8.0 * m6 * P3 - 2.0 * s2_2 * (s2 - 12.0 * m2 * Ta)) / (m * s2_5);
        mu_mu_sigma2_sigma2[i] = (12.0 * s2_3 - 4.0 * s2_3 * DL - 26.0 * m4 * s2 * P2 - 4.0 * m6 * P3 - 32.0 * m2 * s2_2 * Ta) / s2_6;
        mu_sigma2_sigma2_sigma2[i] = 2.0 * (s2_3 * (-14.0 * m + 3.0 * y[i]) + 6.0 * m * s2_3 * DL + 9.0 * m5 * s2 * P2 + m7 * P3 + 18.0 * m3 * s2_2 * Ta) / s2_7;
        sigma2_sigma2_sigma2_sigma2[i] = -(m * (-50.0 * m * s2_3 + 24.0 * s2_3 * y[i] + 24.0 * m * s2_3 * DL + 12.0 * m5 * s2 * P2 + m7 * P3 + 36.0 * m3 * s2_2 * Ta)) / s2_8;
    }

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_sigma2") = mu_mu_mu_sigma2,
        Named("mu_mu_sigma2_sigma2") = mu_mu_sigma2_sigma2,
        Named("mu_sigma2_sigma2_sigma2") = mu_sigma2_sigma2_sigma2,
        Named("sigma2_sigma2_sigma2_sigma2") = sigma2_sigma2_sigma2_sigma2
    );
}

// [[Rcpp::export]]
List gamma_deriv4_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma2) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_sigma2(n), mu_mu_sigma2_sigma2(n),
                  mu_sigma2_sigma2_sigma2(n), sigma2_sigma2_sigma2_sigma2(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool s_is_scalar = (sigma2.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s2 = s_is_scalar ? sigma2[0] : sigma2[i];
        double m2 = m * m, m4 = m2 * m2, m6 = m4 * m2, m8 = m4 * m4;
        double s2_2 = s2 * s2, s2_3 = s2_2 * s2, s2_4 = s2_2 * s2_2,
               s2_5 = s2_4 * s2, s2_6 = s2_3 * s2_3, s2_7 = s2_6 * s2, s2_8 = s2_4 * s2_4;
        double a = m2 / s2;
        double Ta = R::trigamma(a), P2 = R::psigamma(a, 2.0), P3 = R::psigamma(a, 3.0);

        mu_mu_mu_mu[i] = -2.0 * (s2_3 + 8.0 * m6 * P3 + 6.0 * m2 * s2 * (4.0 * m2 * P2 + s2 * Ta)) / (m2 * s2_4);
        mu_mu_mu_sigma2[i] = (36.0 * m4 * s2 * P2 + 8.0 * m6 * P3 - 2.0 * s2_2 * (s2 - 12.0 * m2 * Ta)) / (m * s2_5);
        mu_mu_sigma2_sigma2[i] = -2.0 * (-6.0 * s2_3 + 13.0 * m4 * s2 * P2 + 2.0 * m6 * P3 + 16.0 * m2 * s2_2 * Ta) / s2_6;
        mu_sigma2_sigma2_sigma2[i] = 2.0 * m * (-11.0 * s2_3 + m6 * P3 + 9.0 * m2 * s2 * (m2 * P2 + 2.0 * s2 * Ta)) / s2_7;
        sigma2_sigma2_sigma2_sigma2[i] = (26.0 * m2 * s2_3 - m8 * P3 - 12.0 * m4 * s2 * (m2 * P2 + 3.0 * s2 * Ta)) / s2_8;
    }

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_sigma2") = mu_mu_mu_sigma2,
        Named("mu_mu_sigma2_sigma2") = mu_mu_sigma2_sigma2,
        Named("mu_sigma2_sigma2_sigma2") = mu_sigma2_sigma2_sigma2,
        Named("sigma2_sigma2_sigma2_sigma2") = sigma2_sigma2_sigma2_sigma2
    );
}
