#include <Rcpp.h>
using namespace Rcpp;

// Negative binomial with a variance LINEAR in the mean: Var(Y) = mu (1 + theta),
// against the quadratic mu + mu^2/theta of negbin_distrib(). The two are
// different families rather than two parametrisations of one, and the
// difference shows in where the mean sits: here the size is r = mu/theta, so
// mu appears INSIDE the gamma functions, while in NB2 it stays outside them.
//
// With r = mu/theta and p = 1/(1+theta),
//   l = lgamma(y+r) - lgamma(r) - lgamma(y+1) - r log(1+theta)
//       + y log(theta) - y log(1+theta).
// Writing P = dl/dr and Q = dl/dtheta at fixed r,
//   P    = psi(y+r) - psi(r) - log(1+theta),
//   Q    = -r/(1+theta) + y/theta - y/(1+theta),
//   P_r  = psi'(y+r) - psi'(r),
//   P_th = -1/(1+theta)  (which is also Q_r, the mixed second derivative),
//   Q_th = r/(1+theta)^2 - y/theta^2 + y/(1+theta)^2,
// and the chain rule through r = mu/theta gives the rest.

struct NB1parts {
    double P, Q, Pr, Pth, Qth;
    double rm, rt, rmt, rtt;      // derivatives of r in (mu, theta)
};

static inline NB1parts nb1_parts(double y, double mu, double th) {
    NB1parts z;
    double r = mu / th, om = 1.0 + th;
    z.P   = R::digamma(y + r) - R::digamma(r) - std::log(om);
    z.Q   = -r / om + y / th - y / om;
    z.Pr  = R::trigamma(y + r) - R::trigamma(r);
    z.Pth = -1.0 / om;
    z.Qth = r / (om * om) - y / (th * th) + y / (om * om);
    z.rm  = 1.0 / th;
    z.rt  = -mu / (th * th);
    z.rmt = -1.0 / (th * th);
    z.rtt = 2.0 * mu / (th * th * th);
    return z;
}

// E[psi'(Y + r)] under the family itself, by summing the mass to a far-tail
// quantile. The same device the NB2 kernel uses: there is no closed form, and
// a series against an exact mass is better than a quadrature.
static double nb1_E_trigamma(double mu, double th) {
    double r = mu / th, prob = 1.0 / (1.0 + th);
    double kq = R::qnbinom(1.0 - 1e-12, r, prob, 1, 0);
    int kmax = (int) std::max(100.0, kq) + 1;
    double s = 0.0, cum = 0.0;
    for (int k = 0; k <= kmax; ++k) {
        double pk = R::dnbinom(k, r, prob, 0);
        s += R::trigamma(k + r) * pk;
        cum += pk;
    }
    return (cum > 0) ? s / cum : R::trigamma(r);
}

// [[Rcpp::export]]
NumericVector negbin1_logpmf_cpp(NumericVector y, NumericVector mu,
                                 NumericVector theta) {
    int n = y.size();
    NumericVector out(n);
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);
    for (int i = 0; i < n; i++) {
        double m = mu_s ? mu[0] : mu[i];
        double t = th_s ? theta[0] : theta[i];
        out[i] = R::dnbinom(y[i], m / t, 1.0 / (1.0 + t), 1);
    }
    return out;
}

// [[Rcpp::export]]
List negbin1_gradient_cpp(NumericVector y, NumericVector mu, NumericVector theta) {
    int n = y.size();
    NumericVector g_mu(n), g_th(n);
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);
    for (int i = 0; i < n; i++) {
        double m = mu_s ? mu[0] : mu[i];
        double t = th_s ? theta[0] : theta[i];
        NB1parts z = nb1_parts(y[i], m, t);
        g_mu[i] = z.P * z.rm;
        g_th[i] = z.P * z.rt + z.Q;
    }
    return List::create(Named("mu") = g_mu, Named("theta") = g_th);
}

// [[Rcpp::export]]
List negbin1_hessian_cpp(NumericVector y, NumericVector mu, NumericVector theta) {
    int n = y.size();
    NumericVector h_mm(n), h_mt(n), h_tt(n);
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);
    for (int i = 0; i < n; i++) {
        double m = mu_s ? mu[0] : mu[i];
        double t = th_s ? theta[0] : theta[i];
        NB1parts z = nb1_parts(y[i], m, t);
        h_mm[i] = z.Pr * z.rm * z.rm;
        h_mt[i] = (z.Pr * z.rt + z.Pth) * z.rm + z.P * z.rmt;
        h_tt[i] = z.Pr * z.rt * z.rt + 2.0 * z.Pth * z.rt + z.P * z.rtt + z.Qth;
    }
    return List::create(Named("mu_mu") = h_mm,
                        Named("mu_theta") = h_mt,
                        Named("theta_theta") = h_tt);
}

// The expected Hessian. E[P] = 0 by the first Bartlett identity -- the score
// in mu is P/theta, so its mean vanishing means P's does -- which removes
// every term carrying P and leaves E[psi'(Y+r)] as the only quantity without
// a closed form.
// [[Rcpp::export]]
List negbin1_expected_hessian_cpp(NumericVector y, NumericVector mu,
                                  NumericVector theta) {
    int n = y.size();
    NumericVector h_mm(n), h_mt(n), h_tt(n);
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);
    double last_m = R_NegInf, last_t = R_NegInf, emm = 0, emt = 0, ett = 0;

    for (int i = 0; i < n; i++) {
        double m = mu_s ? mu[0] : mu[i];
        double t = th_s ? theta[0] : theta[i];
        if (m != last_m || t != last_t) {
            double r = m / t, om = 1.0 + t;
            double EPr = nb1_E_trigamma(m, t) - R::trigamma(r);
            double rm = 1.0 / t, rt = -m / (t * t), rtt = 2.0 * m / (t * t * t);
            double Pth = -1.0 / om;
            double EQth = r / (om * om) - m / (t * t) + m / (om * om);
            emm = EPr * rm * rm;
            emt = EPr * rm * rt + Pth * rm;
            ett = EPr * rt * rt + 2.0 * Pth * rt + EQth;
            (void) rtt;   // the term it multiplies has E[P] = 0
            last_m = m; last_t = t;
        }
        h_mm[i] = emm; h_mt[i] = emt; h_tt[i] = ett;
    }
    return List::create(Named("mu_mu") = h_mm,
                        Named("mu_theta") = h_mt,
                        Named("theta_theta") = h_tt);
}
