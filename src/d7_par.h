#ifndef DISTRIBUTIONS7_D7_PAR_H
#define DISTRIBUTIONS7_D7_PAR_H

#include <cstddef>
#include <fenv.h>
#include <RcppParallel.h>

// The one parallel driver of the per-observation kernels. The decomposition
// is over the elements of the OUTPUT: each observation's derivatives are
// computed and written in full by one thread, so no reduction is ever split
// and the result is bit-identical to the sequential loop at any thread
// count. That guarantee is the design rule of the toolkit's parallelism
// (piano_parallel.txt, section 0) and it is what a test holds every kernel
// to: a body that accumulated across observations would not be admissible
// here.
//
// At threads <= 1, or below the kernel's measured threshold, the plain
// sequential loop runs and nothing parallel is entered -- the sequential
// path is the sequential path, not a parallel region with one thread.
// Thresholds are internal and per kernel, measured where the cost of
// opening the region overtakes its gain; they are not user arguments
// (piano_parallel.txt, section 3, on why each of the four knobs a first
// draft had was dropped).
//
// WHAT A BODY MAY CALL. It must not touch the R API: it reads and writes
// preallocated memory through raw pointers, which is what makes it
// admissible inside a worker at all. Two classes of Rmath routine are
// excluded on top of allocation, and they fail differently:
//
//   - anything that can raise a WARNING kills the process. A warning from
//     a worker thread trips R's C-stack check against a foreign stack
//     pointer ("C stack usage ... is too close to the limit", then a
//     segfault) -- observed on four of five CI platforms on 2026-08-19
//     through qnbinom's search reaching pbeta, and reproduced locally on
//     2026-08-21 with lchoose at a non-integer second argument, which
//     warns by design. The p/q family, lchoose/choose and the Bessel
//     functions are all in this class;
//
//   - anything whose value is not thread-invariant breaks the bit-identity
//     the guarantee above promises. Measured locally on 2026-08-21, one
//     input, 120000 calls: psigamma(x, deriv >= 2) at x = 19 and x = 40
//     returns two distinct values 1.3 and 0.8 ulp apart, bessel_k 1.1 ulp,
//     pgamma and pbeta on the log scale 0.5 to 0.8 ulp. Each thread is
//     self-consistent, so which value an element gets follows the chunk it
//     lands in and the answer moves between two runs at the SAME count:
//     gamma1's deriv3 at phi = 1/19 gave six distinct results over six
//     identical calls before this was fixed.
//
// The second class is what fesetenv() below is for, and it is why the
// admissible set may no longer be stated as "the digamma family": digamma
// and trigamma are stable and psigamma at higher orders is not, so the set
// is a measured list and not a family name. Stable at every argument
// probed: digamma, trigamma, lgammafn, lbeta, pnorm, dbinom, dnbinom_mu.
namespace d7 {

// operator() is NOINLINE: the sequential branch and TBB's invocation are
// two call sites, and an inlinable loop gets a separately optimized copy
// at each, which would make the cross-count comparison a comparison of
// two binaries. With the loop out of line, both branches execute the one
// compiled copy; the per-chunk call this costs is nothing against the
// loop it guards.
#if defined(__GNUC__) || defined(__clang__)
#define D7_NOINLINE __attribute__((noinline))
#elif defined(_MSC_VER)
#define D7_NOINLINE __declspec(noinline)
#else
#define D7_NOINLINE
#endif

// The worker carries the CALLING thread's floating-point environment and
// installs it before running its chunk. The control word alone does not
// explain the divergence -- fnstcw reads the same 0x037F on the main
// thread and on a TBB worker -- but restoring the whole environment
// removes it: with it the parallel branch reproduces the sequential value
// exactly at every argument probed, where without it psigamma(19, 2)
// splits 17778 / 182222 over 120000 calls. It costs one fldenv-and-ldmxcsr
// per CHUNK, not per element, and is unmeasurable against the loop.
//
// The environment is captured in the constructor, which runs on the
// calling thread, and installed at the top of operator(), which the
// sequential branch also goes through -- there it restores that thread's
// own environment and is a no-op, so the two branches stay the same code.
template <typename Body>
struct BodyWorker : public RcppParallel::Worker {
  const Body& body;
  fenv_t env;
  explicit BodyWorker(const Body& b) : body(b) { fegetenv(&env); }
  D7_NOINLINE void operator()(std::size_t begin, std::size_t end) {
    fesetenv(&env);
    for (std::size_t i = begin; i < end; ++i) body(i);
  }
};

// The sequential branch runs THROUGH THE WORKER, over the whole range,
// rather than writing a loop of its own: the two branches then execute the
// same compiled function, which is what makes the bit-identity across
// thread counts a property of the code rather than of the optimizer. An
// inlined sequential loop and a worker's out-of-line call are two bodies
// the compiler is free to optimize apart -- observed on the Windows CI
// runner, 2026-08-19, as last-bit differences in the two negbin kernels
// whose per-element arithmetic wraps R::psigamma calls, where nothing
// about the decomposition had changed.
//
// The count is passed to parallelFor rather than left to the process-level
// setting. RcppParallel's resolveValue() gives an explicit positive value
// precedence over RCPP_PARALLEL_NUM_THREADS, so a fit that sized the pool
// through numericals7::local_threads() is unaffected and a caller that did
// not gets the count it asked for: before this, `threads = 2` reached
// parallelFor as -1 and ran on every core the machine has (measured
// 13.9x against the sequential run on a 24-core machine, the same as
// `threads = 24`).
template <typename Body>
inline void par_for(int n, int threads, int threshold, const Body& body) {
  BodyWorker<Body> w(body);
  if (threads > 1 && n >= threshold) {
    RcppParallel::parallelFor(0, static_cast<std::size_t>(n), w, 1, threads);
  } else {
    w(0, static_cast<std::size_t>(n));
  }
}

// Measured crossovers (this machine, 24 cores, threads = 8), each set a
// factor of about two above the measured break-even, where the asymmetry
// argument of the plan's section 3 says an error either way costs almost
// nothing. Region-opening overhead is ~44 microseconds, so the break-even
// moves with the body's cost per observation:
//   - transcendental bodies (digamma and up, ~1 us/obs): parallel already
//     wins 1.3x at n = 64, 3.2x at 256, 9.5x at 2048;
//   - many-output polynomial bodies (third/fourth orders, ~30 ns/obs):
//     0.77x at 2048, 1.1x at 4096, 1.9x at 16384;
//   - two-output polynomial bodies (~8 ns/obs, bandwidth-bound): 0.89x at
//     16384, 1.16x at 32768, 1.5x at 100000.
//   - and a body cheaper still (~4 ns/obs: the geometric's, whose whole
//     score is two divisions) does not break even until about 150000, so it
//     takes a class of its own rather than a threshold picked per family.
//     Measured on the geometric's gradient: 0.79x at 40000, 1.29x at 200000,
//     1.58x at 1000000.
constexpr int kMinCostly = 128;
constexpr int kMinMid = 8192;
constexpr int kMinCheap = 32768;
constexpr int kMinTiny = 131072;

} // namespace d7

#endif
