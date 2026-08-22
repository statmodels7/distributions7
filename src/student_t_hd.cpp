#include <Rcpp.h>
#include "d7_par.h"
using namespace Rcpp;

// Observed third/fourth-order derivatives of the (location-scale) Student's t
// log-density, transcribed from the Wolfram output. r = mu - y,
// D = nu*sigma^2 + r^2. The expected higher derivatives have no closed form
// (handled by the numerical fallback in R).
//
// ⚠️ NOTHING FORMS D OR A POWER OF IT. Every component divided by D^3, and
// D = nu sigma^2 + r^2 overflows at nu sigma^2 = 1.8e+308 while D^3 does so
// at 5.6e+102 -- so at a nu the family's own log link can produce the whole
// surface came back NaN. Measured before this, over nu with sigma = 1: at
// 1e150 three of the ten third derivatives were non-finite, at 1e300 four,
// and at double.xmax eight; on the link scale ten of ten. 0.31.0 made the
// score and the observed Hessian finite to double.xmax by the same rewrite
// and did not reach orders three and four, which is where statmodels7's
// exact outer gradient reads.
//
// The substitution is the one that removes D: with
//
//   z = r/sigma,  u = z^2/nu,  t = 1/(1+u) = nu sigma^2 / D,  a = 1 + 1/nu
//
// every 1/D^k carries a factor (nu sigma^2)^k that the numerator supplies,
// so what is left is a bounded function of (z, u, t) over a power of sigma.
// NINE OF THE TEN ARE EXACT ALGEBRA, with no series and so no crossover to
// calibrate. The two that carry a real cancellation are written out:
//
//   sigma^3 l_sss = 2nu - 2(1+nu) t^3 (1 - 3u), which is -2 in the limit,
//     and equals  -2 + 2 z^2 a (6 + 3u + u^2) t^3  exactly, from
//     1 - t^3(1 - 3u) = u(6 + 3u + u^2) t^3;
//
//   the rational part of l_nunu nu, -4/nu^2 - 8(1+nu)t^3/nu^3 + 12t^2/nu^2,
//     whose bracket -4(1 + 2t^3 - 3t^2) factorizes as -4(t-1)^2(2t+1) with
//     t - 1 = -ut, giving -4 z^4 t^2 (2t+1)/nu^4 - 8 t^3/nu^3.

// psi''((nu+1)/2) - psi''(nu/2) - 8/nu^3 = 12/nu^4 - 20/nu^6 + 84/nu^8 - ...
//
// From the duplication psi(2z) = [psi(z) + psi(z+1/2)]/2 + log 2 twice
// differentiated, psi''(z+1/2) = 8 psi''(2z) - psi''(z), so at z = nu/2 the
// difference is 8 psi''(nu) - 2 psi''(nu/2) and the two Bernoulli expansions
// cancel term by term at nu^-2, leaving 8/nu^3 as the leading behaviour.
// The 8/nu^3 is subtracted HERE, as t_S() carries its own 2/nu^2, because the
// only consumer pairs it with a term that cancels precisely that.
inline double t_T3rest(double v) {
    const double u = 1.0 / v, u2 = u * u, u4 = u2 * u2;
    if (v >= 100.0) return u4 * (12.0 + u2 * (-20.0 + u2 * 84.0));
    return R::psigamma(0.5 * (v + 1.0), 2.0) - R::psigamma(0.5 * v, 2.0) -
        8.0 * u2 * u;
}

// [[Rcpp::export]]
List student_t_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_sigma(n), mu_mu_nu(n), mu_sigma_sigma(n), mu_sigma_nu(n),
                  mu_nu_nu(n), sigma_sigma_sigma(n), sigma_sigma_nu(n), sigma_nu_nu(n), nu_nu_nu(n);
    bool mu_s = (mu.size() == 1), sig_s = (sigma.size() == 1), nu_s = (nu.size() == 1);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_s ? mu[0] : mu[i];
        double s = sig_s ? sigma[0] : sigma[i];
        double v = nu_s ? nu[0] : nu[i];
        double yi = y[i];
        double s2 = s * s, s3 = s2 * s, s4 = s2 * s2;
        double r = m - yi;
        // the dimensionless quantities: nothing below forms D or a power of it
        double iv = 1.0 / v, iv2 = iv * iv, iv3 = iv2 * iv, iv4 = iv2 * iv2;
        double z = r / s, z2 = z * z, z4 = z2 * z2;
        double a = 1.0 + iv;                 // (1 + nu)/nu
        double u = z2 * iv;                  // r^2/(nu sigma^2)
        double t = 1.0 / (1.0 + u);          // nu sigma^2 / D
        double t2 = t * t, t3 = t2 * t;

        mu_mu_mu[i] = 2.0 * a * iv * (3.0 - u) * r * t3 / s4;
        mu_mu_sigma[i] = 2.0 * a * (1.0 - 3.0 * u) * t3 / s3;
        mu_mu_nu[i] = iv2 * (z4 * iv + 1.0 - 3.0 * z2 * a) * t3 / s2;
        mu_sigma_sigma[i] = -2.0 * a * (3.0 - u) * r * t3 / s4;
        mu_sigma_nu[i] = 2.0 * r * iv2 * (z2 * (iv + 2.0) - 1.0) * t3 / s3;
        mu_nu_nu[i] = -2.0 * r * (1.0 - z2) * t3 * iv3 / s2;
        sigma_sigma_sigma[i] =
            (-2.0 + 2.0 * z2 * a * (6.0 + u * (3.0 + u)) * t3) / s3;
        sigma_sigma_nu[i] =
            -z2 * iv2 * (z4 * iv + z2 * (iv + 5.0) - 3.0) * t3 / s2;
        sigma_nu_nu[i] = 2.0 * z2 * (1.0 - z2) * t3 * iv3 / s;
        // 1 - t^3 through expm1/log1p, so the two agree to the last bit
        // however small u is, and t_T3rest carries the 8/nu^3 the second
        // term cancels
        nu_nu_nu[i] = (-4.0 * z4 * iv4 * t2 * (2.0 * t + 1.0)
            - 8.0 * iv3 * std::expm1(-3.0 * std::log1p(u))
            + t_T3rest(v)) / 8.0;
    });

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu, Named("mu_mu_sigma") = mu_mu_sigma, Named("mu_mu_nu") = mu_mu_nu,
        Named("mu_sigma_sigma") = mu_sigma_sigma, Named("mu_sigma_nu") = mu_sigma_nu, Named("mu_nu_nu") = mu_nu_nu,
        Named("sigma_sigma_sigma") = sigma_sigma_sigma, Named("sigma_sigma_nu") = sigma_sigma_nu,
        Named("sigma_nu_nu") = sigma_nu_nu, Named("nu_nu_nu") = nu_nu_nu
    );
}

// [[Rcpp::export]]
List student_t_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_sigma(n), mu_mu_mu_nu(n), mu_mu_sigma_sigma(n),
                  mu_mu_sigma_nu(n), mu_mu_nu_nu(n), mu_sigma_sigma_sigma(n), mu_sigma_sigma_nu(n),
                  mu_sigma_nu_nu(n), mu_nu_nu_nu(n), sigma_sigma_sigma_sigma(n),
                  sigma_sigma_sigma_nu(n), sigma_sigma_nu_nu(n), sigma_nu_nu_nu(n), nu_nu_nu_nu(n);
    bool mu_s = (mu.size() == 1), sig_s = (sigma.size() == 1), nu_s = (nu.size() == 1);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
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
    });

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
