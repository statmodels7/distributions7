# Default Numerical RNG for Continuous Distributions

Fallback method for continuous distributions that do not implement a
native RNG. Two strategies are available and the method picks between
them automatically:

- **Inverse transform**, `distrib_quantile(distrib, runif(n), theta)`,
  when the distribution provides its own quantile function. This is
  exact and costs one quantile evaluation per draw.

- **Generalized Ratio-of-Uniforms**
  ([`rng_grou`](https://statmodels7.github.io/distributions7/reference/rng_grou.md))
  otherwise. Inverting a purely numerical CDF costs one root-finding
  step per draw, which makes simulation-based tools (`approx = "mc"`,
  [`check_distrib`](https://statmodels7.github.io/distributions7/reference/check_distrib.md))
  impractical; GRoU only evaluates the density, so it is orders of
  magnitude faster.

GRoU requires a bounded, unimodal density; if it cannot build its
bounding rectangle the method warns and reverts to inverse transform
sampling. Vector-valued `theta` is handled by grouping the draws by
distinct parameter values.

## Arguments

- distrib:

  An object inheriting from class `"continuous_distrib"`.

- n:

  Number of observations to generate.

- theta:

  A named list of parameters.

## Value

A numeric vector of random draws.
