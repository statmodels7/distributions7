#ifndef D7_PT_SQRT_H
#define D7_PT_SQRT_H

#include <cmath>
#if defined(__SSE2__)
#include <emmintrin.h>
#endif

// The correctly rounded square root at every optimization level. At -O2 gcc
// turns std::sqrt into the sqrtsd instruction, as R's own sqrt() is compiled;
// at -O0, the level pkgload::load_all() builds with, it calls the C library
// instead, and MinGW's sqrt misrounds near-ties (sqrt(0x1.bce198e3f994dp+2)
// gives 0x1.5179af1dfe448p+1 where the root rounds to ...449p+1), so the
// compiled routes stopped being identical() to the R ones on Windows. The
// SSE2 intrinsic is sqrtsd whatever the flags; elsewhere the library sqrt is
// correctly rounded.

namespace d7 {

inline double sqrt_cr(double x) {
#if defined(__SSE2__)
    return _mm_cvtsd_f64(_mm_sqrt_sd(_mm_set_sd(x), _mm_set_sd(x)));
#else
    return std::sqrt(x);
#endif
}

} // namespace d7

#endif
