#ifndef D7_PT_TRANSFORM_H
#define D7_PT_TRANSFORM_H

#include <Rcpp.h>
#include <cmath>
#include <cfloat>

// The ready-made transformers of transformations.R for the scalar registry:
// the inverse x(y) and log|dx/dy| of each, written operation for operation
// as the R closures evaluate them, so that a transformed family's entries
// agree with the R methods to the last bit. R's ^ is R_pow() (powl() on
// 64-bit Windows) and stats' plogis, qlogis and dlogis are nmath's, so both
// are called here; nmath's logistic functions have no warning but the
// silent ME_DOMAIN. The code of a transformer is its position below, and
// tp holds its parameters in the order transformer_scalar_code() lists them.
//
//    1 log         4 sqrt          7 box_cox (lambda)       10 logit
//    2 exp         5 power (p)     8 yeo_johnson (lambda)   11 expit
//    3 inverse     6 asinh         9 affine (loc, scale)    12 softplus (a)

namespace d7 {

const int kTransformN = 12;

// the number of parameters of each transformer, by code
inline int transform_n_par(int code) {
    switch (code) {
    case 5: case 7: case 8: case 12: return 1;
    case 9: return 2;
    default: return 0;
    }
}

// R's sign()
inline double tr_sign(double y) {
    if (std::isnan(y)) return y;
    return (y > 0) ? 1.0 : ((y < 0) ? -1.0 : 0.0);
}

// R's pmax(x, b) for a scalar b, NaN propagating
inline double tr_pmax(double x, double b) {
    return (std::isnan(x) || x >= b) ? x : b;
}

inline double transform_inv(int code, double y, const double* tp) {
    switch (code) {
    case 1: return std::exp(y);
    case 2: return std::log(y);
    case 3: return 1 / y;
    case 4: return y * y;
    case 5: return tr_sign(y) * R_pow(std::fabs(y), 1 / tp[0]);
    case 6: return std::sinh(y);
    case 7: {
        const double lambda = tp[0];
        double base = lambda * y + 1;
        if (base < 0) base = 0;
        return R_pow(base, 1 / lambda);
    }
    case 8: {
        const double lambda = tp[0], lam2 = 2 - lambda;
        if (y >= 0) {
            if (std::fabs(lambda) < 1e-10) return std::expm1(y);
            return R_pow(lambda * y + 1, 1 / lambda) - 1;
        }
        if (std::fabs(lam2) < 1e-10) return -std::expm1(-y);
        return 1 - R_pow(1 - lam2 * y, 1 / lam2);
    }
    case 9: return (y - tp[0]) / tp[1];
    case 10: return R::plogis(y, 0.0, 1.0, 1, 0);
    case 11: return R::qlogis(y, 0.0, 1.0, 1, 0);
    case 12: {
        const double a = tp[0];
        return tr_pmax(y, 0.0) + std::log1p(std::exp(-std::fabs(a * y))) / a;
    }
    default: return R_NaN;
    }
}

inline double transform_log_jac(int code, double y, const double* tp) {
    switch (code) {
    case 1: return y;
    case 2: return -std::log(y);
    case 3: return -2 * std::log(std::fabs(y));
    case 4: return std::log(2.0) + std::log(y);
    case 5: {
        const double p = tp[0];
        return -std::log(std::fabs(p)) + (1 / p - 1) * std::log(std::fabs(y));
    }
    case 6: {
        const double ay = std::fabs(y);
        return ay + std::log1p(std::exp(-2 * ay)) - std::log(2.0);
    }
    case 7: {
        const double lambda = tp[0];
        const double term = tr_pmax(lambda * y + 1, 1e-16);
        return ((1 - lambda) / lambda) * std::log(term);
    }
    case 8: {
        const double lambda = tp[0], lam2 = 2 - lambda;
        if (y >= 0) {
            if (std::fabs(lambda) < 1e-10) return y;
            return ((1 - lambda) / lambda) * std::log1p(lambda * y);
        }
        if (std::fabs(lam2) < 1e-10) return -y;
        return ((lambda - 1) / lam2) * std::log1p(-lam2 * y);
    }
    case 9: return std::log(1 / std::fabs(tp[1]));
    case 10: return R::dlogis(y, 0.0, 1.0, 1);
    case 11: return -(std::log(y) + std::log1p(-y));
    case 12: return R::plogis(tp[0] * y, 0.0, 1.0, 1, 1);
    default: return R_NaN;
    }
}

// the transformed log-density from the parent's at x(y), as
// distrib_pdf(<TransformedDistrib>) composes it: an infinite parent density
// is capped at the largest double, and a zero parent density wins over an
// infinite Jacobian
inline double transform_logpdf(double lpx, double lj) {
    if (std::isinf(lpx) && lpx > 0) lpx = std::log(DBL_MAX);
    if (std::isinf(lpx) && lpx < 0) return R_NegInf;
    return lpx + lj;
}

} // namespace d7

#endif
