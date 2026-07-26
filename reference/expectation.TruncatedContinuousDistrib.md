# Expectation for Truncated Continuous Distributions

The inherited continuous method integrates the density over \\\[\ell,
u\]\\, which is correct unless the parent carries point masses — as it
does when it is a
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
continuous distribution. Those masses are added explicitly, exactly as
in
[`the untruncated case`](https://statmodels7.github.io/distributions7/reference/expectation.ZeroAdjustedContinuousDistrib.md).

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object.

- f:

  A function `f(y, theta, ...)`.

- theta:

  A named list of the parent's parameters.

- ...:

  Additional arguments passed to `f`.
