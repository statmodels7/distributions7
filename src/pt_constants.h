#ifndef D7_PT_CONSTANTS_H
#define D7_PT_CONSTANTS_H

#include <Rcpp.h>

// Constants shared by several families' components.

namespace d7 {

// the Euler-Mascheroni constant
const double kEulerGamma = 0.57721566490153286061;

// (1 - gamma)^2 + pi^2/6, the expected information of the Gumbel's scale
// and of the Weibull's shape, up to a power of the parameter
const double kGumbelInfo = (1.0 - kEulerGamma) * (1.0 - kEulerGamma) +
                           M_PI * M_PI / 6.0;

} // namespace d7

#endif
