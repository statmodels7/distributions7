#ifndef DISTRIBUTIONS7_D7_PAR_H
#define DISTRIBUTIONS7_D7_PAR_H

#include <cstddef>
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
// The body must not touch the R API: it reads and writes preallocated
// memory through raw pointers, which is what makes it admissible inside a
// worker at all.
namespace d7 {

template <typename Body>
struct BodyWorker : public RcppParallel::Worker {
  const Body& body;
  explicit BodyWorker(const Body& b) : body(b) {}
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) body(i);
  }
};

template <typename Body>
inline void par_for(int n, int threads, int threshold, const Body& body) {
  if (threads > 1 && n >= threshold) {
    BodyWorker<Body> w(body);
    RcppParallel::parallelFor(0, static_cast<std::size_t>(n), w);
  } else {
    for (int i = 0; i < n; ++i) body(static_cast<std::size_t>(i));
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
constexpr int kMinCostly = 128;
constexpr int kMinMid = 8192;
constexpr int kMinCheap = 32768;

} // namespace d7

#endif
