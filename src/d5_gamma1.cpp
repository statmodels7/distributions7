#include <Rcpp.h>
#include <cmath>
#include "d7_par.h"
#include "psi_diff.h"
using namespace Rcpp;

// Fifth derivative of the log-density of the gamma in its mean mu and its
// dispersion phi (gamma1 parametrization). The log-density is
//   l = G(a) + a e(mu) - log y,   G(a) = a log a - lgamma(a) - a,
//   e(mu) = log z - (z - 1),      z = y / mu,
// with a = 1 / phi. The fifth derivative is the sum of
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
List gamma1_deriv5_cpp(NumericVector y, NumericVector mu, NumericVector phi,
              int threads = 1) {
    int n = y.size();
    bool m_s = (mu.size() == 1), p_s = (phi.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *pp = phi.begin();
    std::vector<NumericVector> outs;
    std::vector<double*> po(6);
    for (int c = 0; c < 6; ++c) {
        outs.push_back(NumericVector(n));
        po[c] = outs[c].begin();
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        const double m = m_s ? mp[0] : mp[i];
        const double p = p_s ? pp[0] : pp[i];
        const double y = yp[i];
        double res[6];
      const double w = (y - m) / m;
      const double E0 = d7::psi_Ew2(y / m, w);
      const double E1 = (w) / m;
      const double E2 = (( - 1.0 - (2.0*w))) / (m*m);
      const double E3 = ((4.0 + (6.0*w))) / (m*m*m);
      const double E4 = (( - 18.0 - (24.0*w))) / (m*m*m*m);
      const double E5 = ((96.0 + (120.0*w))) / (m*m*m*m*m);
      const double a = 1.0 / p;
      const double G1 = d7::psi_log_rest(a);
      const double G2 = d7::psi1_rest(a);
      const double G3 = d7::psi2_rest(a);
      const double G4 = d7::psi3_rest(a);
      const double G5 = d5_psi4_rest(a);
      const double t0 = (120.0) / (p*p*p*p*p*p);
      res[0] = (E5) / p;
      res[1] = (-1.0*E4) / (p*p);
      res[2] = (2.0*E3) / (p*p*p);
      res[3] = (-6.0*E2) / (p*p*p*p);
      res[4] = (24.0*E1) / (p*p*p*p*p);
      res[5] = ( - (E0*t0) - (G1*t0) - (G5) / (p*p*p*p*p*p*p*p*p*p) - (240.0*G2) / (p*p*p*p*p*p*p) - (120.0*G3) / (p*p*p*p*p*p*p*p) - (20.0*G4) / (p*p*p*p*p*p*p*p*p));
        for (int c = 0; c < 6; ++c) po[c][i] = res[c];
    });
    return List::create(Named("mu_mu_mu_mu_mu") = outs[0], Named("mu_mu_mu_mu_phi") = outs[1], Named("mu_mu_mu_phi_phi") = outs[2], Named("mu_mu_phi_phi_phi") = outs[3], Named("mu_phi_phi_phi_phi") = outs[4], Named("phi_phi_phi_phi_phi") = outs[5]);
}
