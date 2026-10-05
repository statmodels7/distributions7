#ifndef D7_PT_BESSEL_K_H
#define D7_PT_BESSEL_K_H

// numericals7's n7_bessel_k(x, alpha, expo), R's bessel_k() compiled without
// its warnings and without R's allocator, which the pseudo-Huber families
// read for their normalizing constant. It is resolved once through
// R_GetCCallable on the calling thread (bessel_k_fn(), in pseudohuber.cpp);
// the first call must therefore be made there, before any worker starts.

namespace d7 {

typedef double (*N7BesselK)(double, double, double);

N7BesselK bessel_k_fn();

// the exponentially scaled K_alpha(x), as R::bessel_k(x, alpha, 2) returns it
inline double bessel_k_scaled(double x, double alpha) {
    return bessel_k_fn()(x, alpha, 2.0);
}

} // namespace d7

#endif
