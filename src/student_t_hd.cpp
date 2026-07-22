#include <Rcpp.h>
using namespace Rcpp;

// Observed third/fourth-order derivatives of the (location-scale) Student's t
// log-density, transcribed from the Wolfram output. r = mu - y,
// D = nu*sigma^2 + r^2. The expected higher derivatives have no closed form
// (handled by the numerical fallback in R).

// [[Rcpp::export]]
List student_t_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_sigma(n), mu_mu_nu(n), mu_sigma_sigma(n), mu_sigma_nu(n),
                  mu_nu_nu(n), sigma_sigma_sigma(n), sigma_sigma_nu(n), sigma_nu_nu(n), nu_nu_nu(n);
    bool mu_s = (mu.size() == 1), sig_s = (sigma.size() == 1), nu_s = (nu.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = mu_s ? mu[0] : mu[i];
        double s = sig_s ? sigma[0] : sigma[i];
        double v = nu_s ? nu[0] : nu[i];
        double yi = y[i];
        double s2 = s * s, s3 = s2 * s, s4 = s2 * s2, s6 = s4 * s2;
        double r = m - yi, r2 = r * r;
        double D = v * s2 + r2, D2 = D * D, D3 = D2 * D;
        double v1 = 1.0 + v;

        double NUM_mmn = std::pow(m, 4) + v * s4 - 4.0 * std::pow(m, 3) * yi
            + 6.0 * m * v1 * s2 * yi - 3.0 * v1 * s2 * yi * yi
            - 4.0 * m * std::pow(yi, 3) + std::pow(yi, 4)
            - 3.0 * m * m * (v1 * s2 - 2.0 * yi * yi);

        double NUM_ssn = std::pow(m, 4) - 3.0 * v * s4 - 4.0 * std::pow(m, 3) * yi
            + (1.0 + 5.0 * v) * s2 * yi * yi + std::pow(yi, 4)
            - 2.0 * m * yi * ((1.0 + 5.0 * v) * s2 + 2.0 * yi * yi)
            + m * m * ((1.0 + 5.0 * v) * s2 + 6.0 * yi * yi);

        double NUM_msn = m * m * (1.0 + 2.0 * v) - v * s2 + yi * yi + 2.0 * v * yi * yi
            - 2.0 * m * (yi + 2.0 * v * yi);

        mu_mu_mu[i] = -2.0 * v1 * (r2 - 3.0 * v * s2) * r / D3;
        mu_mu_sigma[i] = 2.0 * v * v1 * s * (v * s2 - 3.0 * r2) / D3;
        mu_mu_nu[i] = NUM_mmn / D3;
        mu_sigma_sigma[i] = 2.0 * v * v1 * (r2 - 3.0 * v * s2) * r / D3;
        mu_sigma_nu[i] = 2.0 * s * r * NUM_msn / D3;
        mu_nu_nu[i] = -2.0 * s2 * r * (s2 - r2) / D3;
        sigma_sigma_sigma[i] = 2.0 * v / s3 - 8.0 * std::pow(v, 3) * v1 * s3 / D3 + 6.0 * v * v * v1 * s / D2;
        sigma_sigma_nu[i] = -(r2 * NUM_ssn) / (s2 * D3);
        sigma_nu_nu[i] = 2.0 * s * r2 * (s2 - r2) / D3;
        nu_nu_nu[i] = (-4.0 / (v * v) - 8.0 * v1 * s6 / D3 + 12.0 * s4 / D2
            - R::psigamma(v / 2.0, 2.0) + R::psigamma((1.0 + v) / 2.0, 2.0)) / 8.0;
    }

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu, Named("mu_mu_sigma") = mu_mu_sigma, Named("mu_mu_nu") = mu_mu_nu,
        Named("mu_sigma_sigma") = mu_sigma_sigma, Named("mu_sigma_nu") = mu_sigma_nu, Named("mu_nu_nu") = mu_nu_nu,
        Named("sigma_sigma_sigma") = sigma_sigma_sigma, Named("sigma_sigma_nu") = sigma_sigma_nu,
        Named("sigma_nu_nu") = sigma_nu_nu, Named("nu_nu_nu") = nu_nu_nu
    );
}

// [[Rcpp::export]]
List student_t_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_sigma(n), mu_mu_mu_nu(n), mu_mu_sigma_sigma(n),
                  mu_mu_sigma_nu(n), mu_mu_nu_nu(n), mu_sigma_sigma_sigma(n), mu_sigma_sigma_nu(n),
                  mu_sigma_nu_nu(n), mu_nu_nu_nu(n), sigma_sigma_sigma_sigma(n),
                  sigma_sigma_sigma_nu(n), sigma_sigma_nu_nu(n), sigma_nu_nu_nu(n), nu_nu_nu_nu(n);
    bool mu_s = (mu.size() == 1), sig_s = (sigma.size() == 1), nu_s = (nu.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = mu_s ? mu[0] : mu[i];
        double s = sig_s ? sigma[0] : sigma[i];
        double v = nu_s ? nu[0] : nu[i];
        double yi = y[i];
        double s2 = s * s, s3 = s2 * s, s4 = s2 * s2, s5 = s4 * s, s6 = s4 * s2, s8 = s4 * s4;
        double r = m - yi, r2 = r * r;
        double D = v * s2 + r2, D2 = D * D, D3 = D2 * D, D4 = D2 * D2;
        double v1 = 1.0 + v;
        double m2 = m * m, m3 = m2 * m, m4 = m2 * m2, y2 = yi * yi, y3 = y2 * yi, y4 = y2 * y2;

        // shared with mu_mu_mu_mu and mu_mu_sigma_sigma
        double NUM_m4 = m4 + v * v * s4 - 4.0 * m3 * yi - 6.0 * v * s2 * y2 + y4
            - 4.0 * m * yi * (-3.0 * v * s2 + y2) + 6.0 * m2 * (-(v * s2) + y2);

        double NUM_mmmn = m4 + 3.0 * v * (2.0 + v) * s4 - 4.0 * m3 * yi
            + 4.0 * m * (3.0 + 4.0 * v) * s2 * yi - 2.0 * (3.0 + 4.0 * v) * s2 * y2
            - 4.0 * m * y3 + y4 + m2 * (-2.0 * (3.0 + 4.0 * v) * s2 + 6.0 * y2);

        double NUM_ssnn = m4 + 3.0 * v * s4 - 4.0 * m3 * yi + 2.0 * m * (3.0 + 5.0 * v) * s2 * yi
            - (3.0 + 5.0 * v) * s2 * y2 - 4.0 * m * y3 + y4
            + m2 * (-((3.0 + 5.0 * v) * s2) + 6.0 * y2);

        double NUM_ssnu = m4 * (1.0 + 2.0 * v) + 3.0 * v * v * s4 + 4.0 * m * v * (4.0 + 5.0 * v) * s2 * yi
            - 2.0 * v * (4.0 + 5.0 * v) * s2 * y2 - 4.0 * m * (1.0 + 2.0 * v) * y3 + (1.0 + 2.0 * v) * y4
            - 4.0 * m3 * (yi + 2.0 * v * yi) + 2.0 * m2 * (-(v * (4.0 + 5.0 * v) * s2) + 3.0 * (1.0 + 2.0 * v) * y2);

        double NUM_snn = m4 + v * s4 - 4.0 * m3 * yi + 4.0 * m * v1 * s2 * yi - 2.0 * v1 * s2 * y2
            - 4.0 * m * y3 + y4 - 2.0 * m2 * (v1 * s2 - 3.0 * y2);

        mu_mu_mu_mu[i] = 6.0 * v1 * NUM_m4 / D4;
        mu_mu_mu_sigma[i] = 24.0 * v * v1 * s * (r2 - v * s2) * r / D4;
        mu_mu_mu_nu[i] = -2.0 * r * NUM_mmmn / D4;
        mu_mu_sigma_sigma[i] = -6.0 * v * v1 * NUM_m4 / D4;
        mu_mu_sigma_nu[i] = 2.0 * (-12.0 * v * v * v1 * s5 + 2.0 * v * (7.0 + 9.0 * v) * s3 * D
            - 3.0 * (s + 2.0 * v * s) * D2) / D4;
        mu_mu_nu_nu[i] = 2.0 * (-6.0 * v * v1 * s6 + (5.0 + 9.0 * v) * s4 * D - 3.0 * s2 * D2) / D4;
        mu_sigma_sigma_sigma[i] = -24.0 * v * v * v1 * s * (r2 - v * s2) * r / D4;
        mu_sigma_sigma_nu[i] = 2.0 * r * NUM_ssnu / D4;
        mu_sigma_nu_nu[i] = 4.0 * s * r * NUM_snn / D4;
        mu_nu_nu_nu[i] = 6.0 * s4 * r * (s2 - r2) / D4;
        sigma_sigma_sigma_sigma[i] = 6.0 * v * (-1.0 / s4 + 8.0 * std::pow(v, 3) * v1 * s4 / D4
            - 8.0 * v * v * v1 * s2 / D3 + v * v1 / D2);
        sigma_sigma_sigma_nu[i] = 2.0 / s3 + 24.0 * std::pow(v, 3) * v1 * s5 / D4
            - 4.0 * v * v * (9.0 + 11.0 * v) * s3 / D3 + 6.0 * v * (2.0 + 3.0 * v) * s / D2;
        sigma_sigma_nu_nu[i] = -2.0 * r2 * NUM_ssnn / D4;
        sigma_nu_nu_nu[i] = -6.0 * s3 * r2 * (s2 - r2) / D4;
        nu_nu_nu_nu[i] = 1.0 / (v * v * v) + 3.0 * v1 * s8 / D4 - 4.0 * s6 / D3
            - R::psigamma(v / 2.0, 3.0) / 16.0 + R::psigamma((1.0 + v) / 2.0, 3.0) / 16.0;
    }

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu, Named("mu_mu_mu_sigma") = mu_mu_mu_sigma,
        Named("mu_mu_mu_nu") = mu_mu_mu_nu, Named("mu_mu_sigma_sigma") = mu_mu_sigma_sigma,
        Named("mu_mu_sigma_nu") = mu_mu_sigma_nu, Named("mu_mu_nu_nu") = mu_mu_nu_nu,
        Named("mu_sigma_sigma_sigma") = mu_sigma_sigma_sigma, Named("mu_sigma_sigma_nu") = mu_sigma_sigma_nu,
        Named("mu_sigma_nu_nu") = mu_sigma_nu_nu, Named("mu_nu_nu_nu") = mu_nu_nu_nu,
        Named("sigma_sigma_sigma_sigma") = sigma_sigma_sigma_sigma,
        Named("sigma_sigma_sigma_nu") = sigma_sigma_sigma_nu, Named("sigma_sigma_nu_nu") = sigma_sigma_nu_nu,
        Named("sigma_nu_nu_nu") = sigma_nu_nu_nu, Named("nu_nu_nu_nu") = nu_nu_nu_nu
    );
}
