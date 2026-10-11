#ifndef D7_LS_ORDER_LIST_H
#define D7_LS_ORDER_LIST_H

#include <Rcpp.h>
#include <string>
#include <vector>
#include "d7_par.h"
#include "pt_ls_orders.h"

// The vector kernel of one observed order of a location-scale family in
// (mu, sigma): make(z) returns the family's term functor at the standardized
// residual z, and the components come back in the order of deriv_names(),
// mu_mu..._mu first and sigma_sigma..._sigma last.

namespace d7 {

template <class MakeTerm>
inline Rcpp::List ls_order_list(int order, Rcpp::NumericVector y,
                                Rcpp::NumericVector mu,
                                Rcpp::NumericVector sigma, int threads,
                                MakeTerm make) {
    const int n = y.size();
    std::vector<Rcpp::NumericVector> cols;
    std::vector<double*> p;
    for (int k = 0; k <= order; k++) cols.push_back(Rcpp::NumericVector(n));
    for (int k = 0; k <= order; k++) p.push_back(cols[k].begin());
    const bool mu_is_scalar = (mu.size() == 1);
    const bool sigma_is_scalar = (sigma.size() == 1);

    par_for(n, threads, kMinMid, [&](std::size_t i) {
        const double m = mu_is_scalar ? mu[0] : mu[i];
        const double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double out[6];
        ls_order(order, s, make((y[i] - m) / s), out);
        for (int k = 0; k <= order; k++) p[k][i] = out[k];
    });

    Rcpp::List res(order + 1);
    Rcpp::CharacterVector nm(order + 1);
    for (int k = 0; k <= order; k++) {
        std::string s;
        for (int j = 0; j < order; j++) {
            if (j) s += "_";
            s += (j < order - k) ? "mu" : "sigma";
        }
        res[k] = cols[k];
        nm[k] = s;
    }
    res.attr("names") = nm;
    return res;
}

} // namespace d7

#endif
