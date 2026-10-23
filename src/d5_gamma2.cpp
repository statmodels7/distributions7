#include <Rcpp.h>
#include <cmath>
#include "d7_par.h"
#include "psi_diff.h"
using namespace Rcpp;

// Fifth derivative of the log-density of the gamma in its mean mu and its
// variance sigma2 (gamma2 parametrization). The log-density is
//   l = G(a) + a e(mu) - log y,   G(a) = a log a - lgamma(a) - a,
//   e(mu) = log z - (z - 1),      z = y / mu,
// with a = mu^2 / sigma2. The fifth derivative is the sum of
// G^(r)(a) times the products of the derivatives of a, and of e^(k)(mu) times
// the derivatives of a, written out by the chain rule for the order at hand
// (gamma_gen.py generated the expressions and factored the common
// subexpressions). The shape rests G^(r)(a), r = 1..4, are the psi_diff.h
// forms; G^(5)(a) = -6/a^4 - psi^(4)(a) is d5_psi4_rest below.
// The e-part follows gamma1.cpp: E_0 = [log z - (z - 1)] through psi_Ew2, and
// E_k, k >= 1, from e^(k)(mu) = (-1)^k (k-1)! mu^-k ((1 - k) - k w) with
// w = (y - mu)/mu, the closed form the fourth-order kernels use.

namespace {

// G^(5)(a) = -6/a^4 - psi^(4)(a) = 12/a^5 + 10/a^6 - 7/a^8 + 12/a^10 - 33/a^12 + ...
// For a >= 50 the asymptotic series, whose next term is 130/a^14, is used; the
// truncation is then below 1e-14 relative. Below 50 the direct difference is
// used, where the leading term has not yet dominated.
inline double d5_psi4_rest(double a) {
    if (a >= 50.0) {
        const double u = 1.0 / a, u2 = u * u;
        return u * u * u * u * u *
               (12.0 + u * (10.0 + u2 * (-7.0 + u2 * (12.0 + u2 * (-33.0)))));
    }
    return -6.0 / (a * a * a * a) - R::psigamma(a, 4);
}

} // namespace

// [[Rcpp::export]]
List gamma2_deriv5_cpp(NumericVector y, NumericVector mu, NumericVector sigma2,
              int threads = 1) {
    int n = y.size();
    bool m_s = (mu.size() == 1), p_s = (sigma2.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *pp = sigma2.begin();
    std::vector<NumericVector> outs;
    std::vector<double*> po(6);
    for (int c = 0; c < 6; ++c) {
        outs.push_back(NumericVector(n));
        po[c] = outs[c].begin();
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        const double m = m_s ? mp[0] : mp[i];
        const double v = p_s ? pp[0] : pp[i];
        const double y = yp[i];
        double res[6];
      const double w = (y - m) / m;
      const double E0 = d7::psi_Ew2(y / m, w);
      const double E1 = (w) / m;
      const double E2 = (( - 1.0 - (2.0*w))) / (m*m);
      const double E3 = ((4.0 + (6.0*w))) / (m*m*m);
      const double E4 = (( - 18.0 - (24.0*w))) / (m*m*m*m);
      const double E5 = ((96.0 + (120.0*w))) / (m*m*m*m*m);
      const double a = m * m / v;
      const double G1 = d7::psi_log_rest(a);
      const double G2 = d7::psi1_rest(a);
      const double G3 = d7::psi2_rest(a);
      const double G4 = d7::psi3_rest(a);
      const double G5 = d5_psi4_rest(a);
      const double t0 = (1.0 / (v));
      const double t1 = (1.0 / (v*v*v));
      const double t2 = (t1*m);
      const double t3 = (120.0*G3);
      const double t4 = (m*m);
      const double t5 = (m*m*m);
      const double t6 = (1.0 / (v*v*v*v));
      const double t7 = (m*m*m*m*m);
      const double t8 = (1.0 / (v*v*v*v*v));
      const double t9 = (1.0 / (v*v));
      const double t10 = (12.0*E2);
      const double t11 = (t4*t6);
      const double t12 = (m*m*m*m);
      const double t13 = (m*m*m*m*m*m);
      const double t14 = (1.0 / (v*v*v*v*v*v));
      const double t15 = (t6*m);
      const double t16 = (m*m*m*m*m*m*m);
      const double t17 = (1.0 / (v*v*v*v*v*v*v));
      const double t18 = (12.0*t6);
      const double t19 = (24.0*E1);
      const double t20 = (t4*t8);
      const double t21 = (m*m*m*m*m*m*m*m);
      const double t22 = (1.0 / (v*v*v*v*v*v*v*v));
      const double t23 = (48.0*t8*m);
      const double t24 = (1.0 / (v*v*v*v*v*v*v*v*v));
      const double t25 = (120.0*t14*t4);
      res[0] = ((t2*t3) + (20.0*E3*t0) + (E5*t0*t4) + (10.0*E4*t0*m) + (32.0*G5*t7*t8) + (160.0*G4*t5*t6));
      res[1] = ( - (t10*t9) - (156.0*G3*t11) - (24.0*G2*t1) - (E4*t4*t9) - (112.0*G4*t12*t8) - (16.0*G5*t13*t14) - (8.0*E3*t9*m));
      res[2] = ((t10*t2) + (12.0*E1*t1) + (72.0*G2*t15) + (2.0*E3*t1*t4) + (8.0*G5*t16*t17) + (76.0*G4*t14*t7) + (168.0*G3*t5*t8));
      res[3] = ( - (E0*t18) - (G1*t18) - (t15*t19) - (132.0*G2*t20) - (6.0*E2*t11) - (162.0*G3*t12*t14) - (50.0*G4*t13*t17) - (4.0*G5*t21*t22));
      res[4] = ((E0*t23) + (G1*t23) + (t19*t20) + (2.0*G5*t24*(m*m*m*m*m*m*m*m*m)) + (32.0*G4*t16*t22) + (144.0*G3*t17*t7) + (192.0*G2*t14*t5));
      res[5] = ( - (E0*t25) - (G1*t25) - (G5*(m*m*m*m*m*m*m*m*m*m)) / (v*v*v*v*v*v*v*v*v*v) - (t13*t22*t3) - (240.0*G2*t12*t17) - (20.0*G4*t21*t24));
        for (int c = 0; c < 6; ++c) po[c][i] = res[c];
    });
    return List::create(Named("mu_mu_mu_mu_mu") = outs[0], Named("mu_mu_mu_mu_sigma2") = outs[1], Named("mu_mu_mu_sigma2_sigma2") = outs[2], Named("mu_mu_sigma2_sigma2_sigma2") = outs[3], Named("mu_sigma2_sigma2_sigma2_sigma2") = outs[4], Named("sigma2_sigma2_sigma2_sigma2_sigma2") = outs[5]);
}
