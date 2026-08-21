#include <Rcpp.h>
#include "d7_par.h"
using namespace Rcpp;

// THE DEGREES OF FREEDOM AT LARGE nu.
//
// Three quantities of this family are written as differences that cancel to
// leading order in nu, and each loses every digit it has well before the
// values a fit visits.  Measured on the direct forms:
//
//   A(nu)          rel error 4.1e-08 at nu = 1e4, 0.94 at nu = 1e8
//   T(nu)          rel error 3.9e-11 at nu = 1e5, 9.0e-08 at nu = 1e8
//   E[l_nu_nu]     LOSES ITS SIGN from nu = 3.2e5, reading +2.2e-23 where the
//                  value is -3.5e-24.  On the link scale the factor nu^2 then
//                  makes the information read -500 = -n/2 at nu = 2.7e43, and
//                  a scoring step on a NEGATIVE information walks the fit out
//                  to nu = Inf, where the whole score is NaN.  That is what
//                  made `statmod()` unable to fit this family with nu free.
//
// Each therefore takes an asymptotic branch above a MEASURED crossover, the
// route the generalized Pareto's Lambda(u) already takes here.  The series
// come from the duplication psi(2z) = [psi(z) + psi(z + 1/2)]/2 + log 2 and
// from its derivative, which turn a difference of two psi at nearby arguments
// into a combination at nu and nu/2 whose Bernoulli expansions subtract term
// by term.
//
// The crossovers are where the two routes agree best, measured: A at 200
// (4.3e-10), S at 100 (2.1e-11), E at 1000 (7.0e-08), D at 1e-3 (1.2e-12).
//
// The nu^-6 coefficient of A was WRONG in the first draft (1/8 rather than
// 1/2) and the measurement is what said so: at nu = 100 the series missed the
// direct form by 3.75e-13, which is exactly the difference between the two.

// A(nu) = psi((nu+1)/2) - psi(nu/2) - 1/nu = 1/(2 nu^2) - 1/(4 nu^4) + ...
inline double t_A(double v) {
    if (v >= 200.0) {
        const double u = 1.0 / v, u2 = u * u;
        return u2 * (0.5 + u2 * (-0.25 + u2 * 0.5));
    }
    return R::digamma(0.5 * (v + 1.0)) - R::digamma(0.5 * v) - 1.0 / v;
}

// S(nu) = psi'((nu+1)/2) - psi'(nu/2) + 2/nu^2 = -2/nu^3 + 2/nu^5 - 6/nu^7
// The +2/nu^2 belongs HERE because the only consumer of T pairs it with a
// term that cancels precisely that: forming T alone and cancelling afterwards
// would throw the digits away before the caller sees them.
inline double t_S(double v) {
    const double u = 1.0 / v, u2 = u * u;
    if (v >= 100.0) return u2 * u * (-2.0 + u2 * (2.0 - u2 * 6.0));
    return R::trigamma(0.5 * (v + 1.0)) - R::trigamma(0.5 * v) + 2.0 * u2;
}

// E[l_nu_nu] = [psi'((nu+1)/2) - psi'(nu/2)]/4 + (nu+5)/(2 nu (nu+1)(nu+3))
//            = -7/(2 nu^4) + 13/nu^5 - 79/(2 nu^6) + 119/nu^7 - ...
inline double t_Enunu(double v) {
    if (v >= 1000.0) {
        const double u = 1.0 / v, u2 = u * u, u4 = u2 * u2;
        return u4 * (-3.5 + u * (13.0 + u * (-39.5 + u * 119.0)));
    }
    return 0.25 * (R::trigamma(0.5 * (v + 1.0)) - R::trigamma(0.5 * v)) +
        (v + 5.0) / (2.0 * v * (v + 1.0) * (v + 3.0));
}

// D(u) = u/(1+u) - log1p(u) = -u^2/2 + 2u^3/3 - 3u^4/4 + 4u^5/5 - ...
inline double t_D(double u) {
    if (u < 1e-3) {
        return u * u * (-0.5 + u * (2.0 / 3.0 + u * (-0.75 + u * 0.8)));
    }
    return u / (1.0 + u) - std::log1p(u);
}

// [[Rcpp::export]]
List student_t_gradient_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu,
                        int threads = 1) {
    int n = y.size();
    NumericVector grad_mu(n), grad_sigma(n), grad_nu(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);
    bool nu_is_scalar = (nu.size() == 1);
    bool all_scalar = mu_is_scalar && sigma_is_scalar && nu_is_scalar;

    double m0 = 0, s0 = 0, v0 = 0, s20 = 0, A0 = 0;

    if (all_scalar) {
        m0 = mu[0]; s0 = sigma[0]; v0 = nu[0];
        s20 = s0 * s0;
        A0 = t_A(v0);
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        // LOCAL to the region: a scalar hoisted out of the loop and written
        // inside it is shared once the iterations are split
        double m = m0, s = s0, v = v0, s2 = s20, A = A0;
        if (!all_scalar) {
            m = mu_is_scalar ? mu[0] : mu[i];
            s = sigma_is_scalar ? sigma[0] : sigma[i];
            v = nu_is_scalar ? nu[0] : nu[i];
            s2 = s * s;
            A = t_A(v);
        }
        
        // IN THE RATIO, so that nothing that grows with nu is divided by
        // something else that grows with nu.  Written out, this kernel formed
        // `nu sigma^2`, which OVERFLOWS at the nu the log link can produce
        // (1.8e308 clamped), and then `(nu+1) res / den` is Inf/Inf = NaN --
        // which is how a fit that had legitimately run nu to its boundary
        // came back with a score of NaN and no criterion.  `z2` is formed
        // BEFORE dividing by nu for the same reason.
        double res = y[i] - m;
        double res2 = res * res;
        double z2 = res2 / s2;
        double iv = 1.0 / v;
        double u = z2 * iv;
        double w = 1.0 + u;
        
        grad_mu[i] = (1.0 + iv) * res / (s2 * w);
        grad_sigma[i] = (res2 - s2) / (s * s2 * w);
        // grad_nu = [A(nu) + u/(nu(1+u)) + D(u)] / 2
        // The two terms the direct form writes out, ((nu+1)res^2)/(nu den) and
        // log(1 + res^2/(nu sigma^2)), are both u + O(u^2) and cancel: what
        // survives is D(u) = u/(1+u) - log1p(u), which is -u^2/2 for small u.
        grad_nu[i] = 0.5 * (A + u * iv / w + t_D(u));
    });
    
    return List::create(Named("mu") = grad_mu, Named("sigma") = grad_sigma, Named("nu") = grad_nu);
}

// [[Rcpp::export]]
List student_t_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu,
                        int threads = 1) {
    int n = y.size();
    NumericVector hess_mu_mu(n), hess_sigma_sigma(n), hess_nu_nu(n);
    NumericVector hess_mu_sigma(n), hess_mu_nu(n), hess_sigma_nu(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);
    bool nu_is_scalar = (nu.size() == 1);
    bool all_scalar = mu_is_scalar && sigma_is_scalar && nu_is_scalar;

    double m0 = 0, s0 = 0, v0 = 0, s20 = 0, S0 = 0;

    if (all_scalar) {
        m0 = mu[0]; s0 = sigma[0]; v0 = nu[0];
        s20 = s0 * s0;
        S0 = t_S(v0);
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        // LOCAL to the region, as above
        double m = m0, s = s0, v = v0, s2 = s20, S = S0;
        if (!all_scalar) {
            m = mu_is_scalar ? mu[0] : mu[i];
            s = sigma_is_scalar ? sigma[0] : sigma[i];
            v = nu_is_scalar ? nu[0] : nu[i];
            s2 = s * s;
            S = t_S(v);
        }
        
        // IN THE RATIO throughout -- see the gradient above for why.
        double res = y[i] - m;
        double res2 = res * res;
        double res4 = res2 * res2;
        double z2 = res2 / s2;
        double iv = 1.0 / v;
        double iv2 = iv * iv;
        double u = z2 * iv;
        double w = 1.0 + u;
        double w2 = w * w;
        
        hess_mu_mu[i] = (1.0 + iv) * (u - 1.0) / (s2 * w2);
        // numerator v(v s^4 - (3v+1) s^2 res^2 - res^4) = v^2 s^4 (1 - (3v+1)u - v u^2),
        // and with v u = z^2 that is 1 - 3 z^2 - u - z^2 u
        hess_sigma_sigma[i] = (1.0 - 3.0 * z2 - u - z2 * u) / (s2 * w2);
        // hess_nu_nu = [S(nu) - (2/nu^2) u(2+u)/(1+u)^2 + 2 z^4/(nu^3 (1+u)^2)]/4
        // S carries T(nu) + 2/nu^2 because the 2/nu^2 of the second term is
        // exactly what cancels T's leading order; the remainder of that term
        // is written as -(2/nu^2)(1/(1+u)^2 - 1), which has no cancellation.
        hess_nu_nu[i] = 0.25 * (S - 2.0 * iv2 * u * (2.0 + u) / w2
                                + 2.0 * z2 * z2 * iv2 * iv / w2);
        
        hess_mu_sigma[i] = -2.0 * (1.0 + iv) * res / (s * s2 * w2);
        hess_mu_nu[i] = res * (res2 - s2) * iv2 / (s2 * s2 * w2);
        hess_sigma_nu[i] = (res4 - s2 * res2) * iv2 / (s * s2 * s2 * w2);
    });
    
    return List::create(
        Named("mu_mu") = hess_mu_mu, Named("sigma_sigma") = hess_sigma_sigma, Named("nu_nu") = hess_nu_nu,
        Named("mu_sigma") = hess_mu_sigma, Named("mu_nu") = hess_mu_nu, Named("sigma_nu") = hess_sigma_nu
    );
}

// [[Rcpp::export]]
List student_t_expected_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu,
                        int threads = 1) {
    int n = y.size();
    NumericVector hess_mu_mu(n), hess_sigma_sigma(n), hess_nu_nu(n);
    NumericVector hess_mu_sigma(n), hess_mu_nu(n), hess_sigma_nu(n);
    
    bool sigma_is_scalar = (sigma.size() == 1);
    bool nu_is_scalar = (nu.size() == 1);
    bool both_scalar = sigma_is_scalar && nu_is_scalar;

    double s0 = 0, v0 = 0, s20 = 0, E0 = 0;

    if (both_scalar) {
        s0 = sigma[0]; v0 = nu[0];
        s20 = s0 * s0;
        E0 = t_Enunu(v0);
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        // LOCAL to the region, as above
        double s = s0, v = v0, s2 = s20, E = E0;
        if (!both_scalar) {
            s = sigma_is_scalar ? sigma[0] : sigma[i];
            v = nu_is_scalar ? nu[0] : nu[i];
            s2 = s * s;
            E = t_Enunu(v);
        }
        
        // IN THE RATIO: (v+1)/(v+3) is Inf/Inf at the nu the log link can
        // produce, and reads NaN where the value is 1.
        double iv = 1.0 / v;
        hess_mu_mu[i] = -(1.0 + iv) / (s2 * (1.0 + 3.0 * iv));
        hess_sigma_sigma[i] = -2.0 / (s2 * (1.0 + 3.0 * iv));
        hess_nu_nu[i] = E;
        
        hess_sigma_nu[i] = 2.0 * iv * iv / ((1.0 + iv) * (1.0 + 3.0 * iv) * s);
        // mu_sigma and mu_nu default initialized to 0.0, no need to assign
    });
    
    return List::create(
        Named("mu_mu") = hess_mu_mu, Named("sigma_sigma") = hess_sigma_sigma, Named("nu_nu") = hess_nu_nu,
        Named("mu_sigma") = hess_mu_sigma, Named("mu_nu") = hess_mu_nu, Named("sigma_nu") = hess_sigma_nu
    );
}