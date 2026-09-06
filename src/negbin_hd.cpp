#include <Rcpp.h>
#include "d7_par.h"
using namespace Rcpp;

// Third/fourth-order derivatives of the Negative Binomial (NB2) log-mass,
// derived from the gradient/Hessian by repeated differentiation and validated
// numerically. s = theta + mu.
//
// The expected pure-theta derivatives need an expectation of a polygamma over
// the support, summed with the same pmf recurrence as the expected Hessian's
// nb_E_trigamma_diff, carried in log scale until p_k is representable and
// stopped when the accumulated mass reaches 1 - 1e-12. No quantile or mass
// function is called: this helper runs inside d7::par_for workers, and
// qnbinom's search reaches pbeta, whose warning path calls into the R API and
// killed the process from a worker thread on four of the five CI platforms.
//
// Both the OBSERVED and the EXPECTED pure-theta derivatives are written so
// that the cancellation is performed symbolically; the two devices are below.
// THE PURE-theta DERIVATIVES AT ORDERS THREE AND FOUR, OBSERVED AND EXPECTED.
//
// Each is a polygamma difference paired with closed terms that cancel it to
// leading order as theta runs away, exactly as the score and the hessian do
// one file over.  Measured against the forms below, which perform each
// cancellation symbolically: the spelling replaced is 4.15e-02 out on the
// observed third derivative at theta = 1e7 and 3.09e-02 on the fourth, and the
// composed expected ones LOSE THEIR SIGN there -- -3.23e-32 where the value is
// +4.80e-34 at order three, +9.55e-39 where it is -2.88e-40 at order four.
//
// THE OBSERVED ONES.  The shift is y, a count, so
//   psi''(y+th) - psi''(th)  =  2 sum_{j<y} 1/(th+j)^3
//   psi'''(y+th) - psi'''(th) = -6 sum_{j<y} 1/(th+j)^4
// and the term in (y - mu)/s^k that accompanies each is itself a sum over the
// same range plus a constant, since y = sum_{j<y} 1.  Merging them and
// factorizing, with t = th + j and s = th + mu,
//
//   l_ttt  = sum_{j<y} 2(mu-j)(s^2 + s t + t^2)/(t^3 s^3)
//            - mu^2 (3th + mu)/(th^2 s^3)
//   l_tttt = -sum_{j<y} 6(mu-j)(t+s)(t^2+s^2)/(t^4 s^4)
//            + 2 mu^2 (6th^2 + 4 th mu + mu^2)/(th^3 s^4)
//
// every factor an exact polynomial.  The sum costs y terms, so it is taken
// only where the direct form loses digits, which is where theta is large
// against y -- and there the sum is SHORT, y being below th/kObsCut.  Measured,
// the direct form is within 3.9e-10 at th/y <= 100 and reaches 5e-02 at 1e6,
// so the threshold is 100.  On a real design of 57600 counts summing to
// 6.24e6 the whole sum is 6.24e6 terms at the worst theta, which is a few
// milliseconds.
//
// THE EXPECTED ONES.  The closed term is mu times a constant, and
// E[sum_{j<Y} c] = c E[Y] = c mu, so it is itself an expectation over the same
// support and merges into the summand rather than being added afterwards.
// With K = 3th^2 + 3 th mu + mu^2 the merged terms are
//
//   order 3:  [ th^2 mu (3th + 2mu) - (2th+mu) j (3th^2 + 3 j th + j^2) ]
//             / (t^3 th^2 s^2)
//   order 4:  [ -2 th^3 mu (6th^2 + 8 th mu + 3mu^2)
//               + 2 K j (4th^3 + 6 j th^2 + 4 j^2 th + j^3) ] / (t^4 th^3 s^3)
//
// both exact against the naive term to 3e-16 at every probe tried.
// ⚠️ Neither device removes the cancellation, each moves it from order theta
// to order theta/mu, which is the same trade nb_E_ltt records one file over.

static const double kObsCut = 100.0;   // use the sum when th > kObsCut * y

// l_theta^3 and l_theta^4 at one observation, without the cancellation
static inline double nb_d3_theta(double y, double mu, double th) {
    const double s = th + mu, s3 = s * s * s;
    if (!(th > kObsCut * (y > 1.0 ? y : 1.0))) {
        const double s2 = s * s;
        return R::psigamma(y + th, 2.0) - R::psigamma(th, 2.0)
            - mu * (2.0 * th + mu) / (th * th * s2) - 2.0 * (y - mu) / s3;
    }
    double acc = 0.0;
    for (double j = 0.0; j < y; j += 1.0) {
        const double t = th + j;
        acc += 2.0 * (mu - j) * (s * s + s * t + t * t) / (t * t * t * s3);
    }
    return acc - mu * mu * (3.0 * th + mu) / (th * th * s3);
}

static inline double nb_d4_theta(double y, double mu, double th) {
    const double s = th + mu, s2 = s * s, s4 = s2 * s2;
    if (!(th > kObsCut * (y > 1.0 ? y : 1.0))) {
        const double s3 = s2 * s;
        return R::psigamma(y + th, 3.0) - R::psigamma(th, 3.0)
            - 2.0 * mu * (th * s - (2.0 * th + mu) * (2.0 * th + mu))
              / (th * th * th * s3)
            + 6.0 * (y - mu) / s4;
    }
    double acc = 0.0;
    for (double j = 0.0; j < y; j += 1.0) {
        const double t = th + j, t2 = t * t;
        acc -= 6.0 * (mu - j) * (t + s) * (t2 + s2) / (t2 * t2 * s4);
    }
    return acc + 2.0 * mu * mu * (6.0 * th * th + 4.0 * th * mu + mu * mu)
                 / (th * th * th * s4);
}

// E[l_theta^3] and E[l_theta^4], each as ONE sum over the support.  The mass
// comes from the pmf recurrence, seeded and switched exactly as
// nb_E_trigamma's is; see the head of negbin.cpp for both rules.
static double nb_E_dtheta(double mu, double theta, int nd) {
    const double ratio = mu / (theta + mu);
    const double lratio = std::log(mu) - std::log(theta + mu);
    const double cap = 100.0 + mu + 20.0 * std::sqrt(mu * (1.0 + mu / theta))
                       + 40.0 * (mu + theta) / theta;
    const int kmax = (int) std::min(cap, 2.0e9);
    const double s = theta + mu, th2 = theta * theta, th3 = th2 * theta;
    const double s2 = s * s, s3 = s2 * s;
    const double A3 = th2 * mu * (3.0 * theta + 2.0 * mu);
    const double B3 = 2.0 * theta + mu;
    const double A4 = -2.0 * th3 * mu
        * (6.0 * th2 + 8.0 * theta * mu + 3.0 * mu * mu);
    const double K4 = 3.0 * th2 + 3.0 * theta * mu + mu * mu;

    double acc = 0.0, cum = 0.0, U = 0.0;
    double lpk = -theta * std::log1p(mu / theta);
    bool logscale = (lpk <= -640.0);
    double pk = std::exp(lpk);
    int k = 0;
    for (; k <= kmax; ++k) {
        acc += U * pk;
        cum += pk;
        const bool last = (cum >= 1.0 - 1e-12 && k >= 100);
        const double j = (double) k, t = theta + j, t2 = t * t;
        if (nd == 2) {
            U += (A3 - B3 * j * (3.0 * th2 + 3.0 * j * theta + j * j))
                 / (t2 * t * th2 * s2);
        } else {
            U += (A4 + 2.0 * K4 * j
                  * (4.0 * th3 + 6.0 * j * th2 + 4.0 * j * j * theta + j * j * j))
                 / (t2 * t2 * th3 * s3);
        }
        if (last) { ++k; break; }
        if (logscale) {
            lpk += lratio + std::log((k + theta) / (k + 1.0));
            pk = std::exp(lpk);
            if (lpk > -640.0) logscale = false;
        } else {
            pk *= (k + theta) / (k + 1.0) * ratio;
        }
    }
    if (cum < 1.0) acc += U * (1.0 - cum);
    return acc;
}

// [[Rcpp::export]]
List negbin_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector theta,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_theta(n), mu_theta_theta(n), theta_theta_theta(n);
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_s ? mu[0] : mu[i];
        double th = th_s ? theta[0] : theta[i];
        double yi = y[i];
        double s = th + m, s2 = s * s, s3 = s2 * s;
        double m2 = m * m, m3 = m2 * m;

        mu_mu_mu[i] = -2.0 * (yi + th) / s3 + 2.0 * yi / m3;
        mu_mu_theta[i] = 1.0 / s2 - 2.0 * (yi + th) / s3;
        mu_theta_theta[i] = -2.0 * (yi - m) / s3;
        theta_theta_theta[i] = nb_d3_theta(yi, m, th);
    });

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_theta") = mu_mu_theta,
        Named("mu_theta_theta") = mu_theta_theta,
        Named("theta_theta_theta") = theta_theta_theta
    );
}

// [[Rcpp::export]]
List negbin_deriv3_expected_cpp(NumericVector y, NumericVector mu, NumericVector theta,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_theta(n), mu_theta_theta(n), theta_theta_theta(n);
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);
    bool both_scalar = mu_s && th_s;

    // the scalar-case constant lives OUT here; the per-iteration copy is
    // LOCAL to the lambda, or two threads would race on it
    double ttt0 = 0;
    if (both_scalar) {
        double m = mu[0], th = theta[0];
        ttt0 = nb_E_dtheta(m, th, 2);
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_s ? mu[0] : mu[i];
        double th = th_s ? theta[0] : theta[i];
        double s = th + m, s2 = s * s;
        double ttt = ttt0;
        if (!both_scalar) {
            ttt = nb_E_dtheta(m, th, 2);
        }
        mu_mu_mu[i] = -2.0 / s2 + 2.0 / (m * m);
        mu_mu_theta[i] = -1.0 / s2;
        mu_theta_theta[i] = 0.0;
        theta_theta_theta[i] = ttt;
    });

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_theta") = mu_mu_theta,
        Named("mu_theta_theta") = mu_theta_theta,
        Named("theta_theta_theta") = theta_theta_theta
    );
}

// [[Rcpp::export]]
List negbin_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector theta,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_theta(n), mu_mu_theta_theta(n),
                  mu_theta_theta_theta(n), theta_theta_theta_theta(n);
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_s ? mu[0] : mu[i];
        double th = th_s ? theta[0] : theta[i];
        double yi = y[i];
        double s = th + m, s2 = s * s, s3 = s2 * s, s4 = s2 * s2;
        double m2 = m * m, m4 = m2 * m2;

        mu_mu_mu_mu[i] = 6.0 * (yi + th) / s4 - 6.0 * yi / m4;
        mu_mu_mu_theta[i] = -2.0 / s3 + 6.0 * (yi + th) / s4;
        mu_mu_theta_theta[i] = -4.0 / s3 + 6.0 * (yi + th) / s4;
        mu_theta_theta_theta[i] = 6.0 * (yi - m) / s4;
        theta_theta_theta_theta[i] = nb_d4_theta(yi, m, th);
    });

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_theta") = mu_mu_mu_theta,
        Named("mu_mu_theta_theta") = mu_mu_theta_theta,
        Named("mu_theta_theta_theta") = mu_theta_theta_theta,
        Named("theta_theta_theta_theta") = theta_theta_theta_theta
    );
}

// [[Rcpp::export]]
List negbin_deriv4_expected_cpp(NumericVector y, NumericVector mu, NumericVector theta,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_theta(n), mu_mu_theta_theta(n),
                  mu_theta_theta_theta(n), theta_theta_theta_theta(n);
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);
    bool both_scalar = mu_s && th_s;

    double tttt0 = 0;
    if (both_scalar) {
        double m = mu[0], th = theta[0];
        tttt0 = nb_E_dtheta(m, th, 3);
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_s ? mu[0] : mu[i];
        double th = th_s ? theta[0] : theta[i];
        double s = th + m, s3 = s * s * s;
        double tttt = tttt0;
        if (!both_scalar) {
            tttt = nb_E_dtheta(m, th, 3);
        }
        mu_mu_mu_mu[i] = 6.0 / s3 - 6.0 / (m * m * m);
        mu_mu_mu_theta[i] = 4.0 / s3;
        mu_mu_theta_theta[i] = 2.0 / s3;
        mu_theta_theta_theta[i] = 0.0;
        theta_theta_theta_theta[i] = tttt;
    });

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_theta") = mu_mu_mu_theta,
        Named("mu_mu_theta_theta") = mu_mu_theta_theta,
        Named("mu_theta_theta_theta") = mu_theta_theta_theta,
        Named("theta_theta_theta_theta") = theta_theta_theta_theta
    );
}
