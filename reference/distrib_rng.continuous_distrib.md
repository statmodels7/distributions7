# Default Random Generation for Continuous Distributions

The fallback for a continuous family that implements no generator of its
own. It picks between two routes by asking whether the family has a real
quantile method:

- **inverse transform** where
  [`has_analytic_quantile()`](https://statmodels7.github.io/distributions7/reference/has_analytic_quantile.md)
  is `TRUE`, which costs one quantile evaluation per draw;

- **generalized ratio-of-uniforms**
  ([`rng_grou()`](https://statmodels7.github.io/distributions7/reference/rng_grou.md))
  otherwise, which evaluates only the density.

The choice is what keeps simulation-based tools usable at all. Inverting
a purely numerical distribution function costs one root-finding step per
draw, and each of those is itself a quadrature, so `approx = "mc"` and
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
would be impractical on a family defined by its density alone. Measured
on such a family, 2000 draws take under a millisecond by this route and
pass a Kolmogorov-Smirnov test against the truth.

## Arguments

- distrib:

  An object inheriting from `continuous_distrib` that registers no
  method of its own.

- n:

  The number of draws, a single non-negative integer.

- theta:

  A named list of parameters, each of length 1 or `n`. A component of
  length `n` gives one draw per parameter setting.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `n` draws.

## Details

The ratio-of-uniforms scheme needs a bounded unimodal density. Where its
bounding rectangle cannot be built the method **warns and reverts** to
inverse transform sampling, so a draw is always returned. A `theta`
varying by observation is handled by grouping the draws by distinct
parameter setting, one bounding rectangle per group.

## See also

[`rng_grou()`](https://statmodels7.github.io/distributions7/reference/rng_grou.md)
for the ratio-of-uniforms sampler and its theory;
[`has_analytic_quantile()`](https://statmodels7.github.io/distributions7/reference/has_analytic_quantile.md),
which chooses between the two routes;
[`distrib_rng.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.discrete_distrib.md),
where the question does not arise;
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.
