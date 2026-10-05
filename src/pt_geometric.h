#ifndef D7_PT_GEOMETRIC_H
#define D7_PT_GEOMETRIC_H

// Geometric on {0, 1, ...} in its mean, l = y log mu - (y + 1) log(1 + mu):
// one function per component and order for the quantities the scalar
// registry reads, called by the vector kernels (geometric.cpp,
// dexpected_kernels.cpp) and by the registry (d7_ccallable.cpp).

namespace d7 {

inline double geometric_score_mu(double y, double m) {
    return (y - m) / (m * (1.0 + m));
}

inline double geometric_hess_mu_mu(double y, double m) {
    double om = 1.0 + m;
    return -(y / (m * m) - (y + 1.0) / (om * om));
}

inline double geometric_expected_mu_mu(double m) {
    return -1.0 / (m * (1.0 + m));
}

// E = -1/q with q = m(1 + m)
inline double geometric_dexpected_mu_mu_mu(double m) {
    double q = m * (1.0 + m), s = 1.0 + 2.0 * m, q2 = q * q;
    return s / q2;
}

inline void geometric_score_curv(int k, double y, const double* th,
                                 double* out) {
    out[0] = geometric_score_mu(y, th[0]);
    out[1] = geometric_hess_mu_mu(y, th[0]);
}

inline void geometric_info_dinfo(int k, double y, const double* th,
                                 double* out) {
    out[0] = geometric_expected_mu_mu(th[0]);
    out[1] = geometric_dexpected_mu_mu_mu(th[0]);
}

} // namespace d7

#endif
