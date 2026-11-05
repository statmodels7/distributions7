// The derivatives of the log density on the link scale, from those on the
// parameter scale and the derivatives of the inverse links, at one
// observation (piano_motore.md, D12). Header only, read through
// `LinkingTo: distributions7`.
//
// Each parameter has its own link, theta_i = h_i(eta_i), so the Jacobian is
// diagonal and the multivariate Faa di Bruno formula factorizes over the
// distinct parameters of a multi-index. Every function below repeats the
// operations of to_link_scale() in R/link_scale.R in the same order, so that
// the two routes agree bit for bit: the accumulator starts at 0, the terms
// are added in the order of link_scale_layout()'s combinations, and each
// coefficient is the product, in parameter order, of 1 and the partial Bell
// polynomials of bell_partial().
//
// Layout: K parameters; g[i] the score in theta_i; H the K x K Hessian in
// theta, column-major and complete; h1[i], h2[i] the first and second
// derivatives of h_i at eta_i. The outputs are complete K x K matrices,
// column-major. One function per order.
#ifndef D7_LINK_SCALE_H
#define D7_LINK_SCALE_H

namespace d7ls {

// order 1: dl/deta_i = (dl/dtheta_i) h_i'
inline void link_grad(int K, const double* g, const double* h1, double* out) {
    for (int i = 0; i < K; ++i) out[i] = g[i] * h1[i];
}

// order 2 of the observed Hessian. On the diagonal the two combinations of
// the layout are the score with B_{2,1} = h'' and the second derivative with
// B_{2,2} = h'^2; off the diagonal the single combination has the
// coefficient (1 * h_i') * h_j'.
inline void link_hess(int K, const double* g, const double* H,
                      const double* h1, const double* h2, double* out) {
    for (int i = 0; i < K; ++i) {
        double acc = 0.0;
        acc = acc + g[i] * (1.0 * h2[i]);
        acc = acc + H[i + i * K] * (1.0 * (h1[i] * h1[i]));
        out[i + i * K] = acc;
        for (int j = i + 1; j < K; ++j) {
            double v = 0.0;
            v = v + H[i + j * K] * ((1.0 * h1[i]) * h1[j]);
            out[i + j * K] = v;
            out[j + i * K] = v;
        }
    }
}

// order 2 of the expected Hessian: the generic passes a zero score
// (distrib_expected_hessian(), `zero1`), so the first combination is
// 0 * h'' and is kept, since it is NaN where h'' is not finite
inline void link_info(int K, const double* E, const double* h1,
                      const double* h2, double* out) {
    for (int i = 0; i < K; ++i) {
        double acc = 0.0;
        acc = acc + 0.0 * (1.0 * h2[i]);
        acc = acc + E[i + i * K] * (1.0 * (h1[i] * h1[i]));
        out[i + i * K] = acc;
        for (int j = i + 1; j < K; ++j) {
            double v = 0.0;
            v = v + E[i + j * K] * ((1.0 * h1[i]) * h1[j]);
            out[i + j * K] = v;
            out[j + i * K] = v;
        }
    }
}

}  // namespace d7ls

#endif
