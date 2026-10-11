#include <Rcpp.h>
#include <string>
#include <vector>
#include "d7_par.h"
#include "pt_cauchy.h"
using namespace Rcpp;

// Third, fourth and fifth observed derivatives of the Cauchy log-density, and
// the derivatives of its distribution function to fourth order, from the
// closed forms of pt_cauchy.h. Components in the order of deriv_names():
// mu_..._mu first, sigma_..._sigma last.

static List cauchy_list(int order, NumericVector y, NumericVector mu,
                        NumericVector sigma, int threads, bool cdf) {
    const int n = y.size();
    std::vector<NumericVector> cols;
    std::vector<double*> p;
    for (int k = 0; k <= order; k++) cols.push_back(NumericVector(n));
    for (int k = 0; k <= order; k++) p.push_back(cols[k].begin());
    const bool mu_is_scalar = (mu.size() == 1);
    const bool sigma_is_scalar = (sigma.size() == 1);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        const double m = mu_is_scalar ? mu[0] : mu[i];
        const double s = sigma_is_scalar ? sigma[0] : sigma[i];
        const double res = y[i] - m;
        for (int k = 0; k <= order; k++) {
            p[k][i] = cdf ? d7::cauchy_cdf_d(order, order - k, res, s)
                          : d7::cauchy_logdens_d(order, order - k, res, s);
        }
    });

    List out(order + 1);
    CharacterVector nm(order + 1);
    for (int k = 0; k <= order; k++) {
        std::string s;
        for (int j = 0; j < order; j++) {
            if (j) s += "_";
            s += (j < order - k) ? "mu" : "sigma";
        }
        out[k] = cols[k];
        nm[k] = s;
    }
    out.attr("names") = nm;
    return out;
}

// [[Rcpp::export]]
List cauchy_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                        int threads = 1) {
    return cauchy_list(3, y, mu, sigma, threads, false);
}

// [[Rcpp::export]]
List cauchy_deriv3_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_sigma(n), mu_sigma_sigma(n), sigma_sigma_sigma(n);
    bool sigma_is_scalar = (sigma.size() == 1);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double s3 = s * s * s;
        mu_mu_mu[i] = 0.0;
        mu_mu_sigma[i] = 1.0 / (2.0 * s3);
        mu_sigma_sigma[i] = 0.0;
        sigma_sigma_sigma[i] = 3.0 / (2.0 * s3);
    });

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_sigma") = mu_mu_sigma,
        Named("mu_sigma_sigma") = mu_sigma_sigma,
        Named("sigma_sigma_sigma") = sigma_sigma_sigma
    );
}

// [[Rcpp::export]]
List cauchy_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                        int threads = 1) {
    return cauchy_list(4, y, mu, sigma, threads, false);
}

// [[Rcpp::export]]
List cauchy_deriv4_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_sigma(n), mu_mu_sigma_sigma(n),
                  mu_sigma_sigma_sigma(n), sigma_sigma_sigma_sigma(n);
    bool sigma_is_scalar = (sigma.size() == 1);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double s4 = s * s * s * s;
        mu_mu_mu_mu[i] = 3.0 / (4.0 * s4);
        mu_mu_mu_sigma[i] = 0.0;
        mu_mu_sigma_sigma[i] = -3.0 / (4.0 * s4);
        mu_sigma_sigma_sigma[i] = 0.0;
        sigma_sigma_sigma_sigma[i] = -21.0 / (4.0 * s4);
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
List cauchy_deriv5_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                        int threads = 1) {
    return cauchy_list(5, y, mu, sigma, threads, false);
}

// the order-k derivatives of the distribution function, lower tail
// [[Rcpp::export]]
List cauchy_cdf_deriv_cpp(NumericVector q, NumericVector mu, NumericVector sigma,
                          int order, int threads = 1) {
    return cauchy_list(order, q, mu, sigma, threads, true);
}
