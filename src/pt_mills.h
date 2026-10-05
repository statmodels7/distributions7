#ifndef D7_PT_MILLS_H
#define D7_PT_MILLS_H

#include <Rcpp.h>
#include <cmath>

namespace d7 {

// The Mills ratio phi(t)/Phi(t), on the log scale so that both factors may
// underflow in the left tail while their ratio stays finite; it is
// numericals7::mills_ratio() written in C. Its derivative is -r (t + r).
inline double mills_ratio(double t) {
    return std::exp(R::dnorm4(t, 0.0, 1.0, 1) - R::pnorm5(t, 0.0, 1.0, 1, 1));
}

} // namespace d7

#endif
