# Default Numerical CDF for Continuous Distributions

The fallback for a continuous family that implements no analytical
distribution function: it integrates `distrib_pdf` numerically. Measured
on a Gamma defined by its density alone, it agrees with
[`stats::pgamma()`](https://rdrr.io/r/stats/GammaDist.html) to between
\\1.2\times10^{-16}\\ and \\2.6\times10^{-15}\\ relative.

## Arguments

- distrib:

  An object inheriting from `continuous_distrib` that registers no
  method of its own.

- q:

  A numeric vector of quantiles. A value outside `distrib@bounds`
  returns 0 or 1 as the side requires.

- theta:

  A named list of parameters, each of length 1 or `length(q)`.

- lower.tail:

  Logical of length 1. `TRUE`, the default, returns \\P(Y \le q)\\;
  `FALSE` returns the upper tail, computed as the complement.

- log.p:

  Logical of length 1. `TRUE` returns the logarithm of the probability.
  The logarithm is taken after the integration, so it does not recover a
  tail that has already underflowed; a family whose far tail matters
  should write its own method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of cumulative probabilities, the length of the recycled
`q` and `theta`. `NA` for a row whose quadrature did not reach its
accuracy target.

## Where the nodes go

An approximate mode is located first by
[`find_pdf_anchor()`](https://statmodels7.github.io/distributions7/reference/find_pdf_anchor.md),
and the integral is taken over whichever side of the mode holds \\q\\,
the other side being reached through the complement. The quadrature
nodes then concentrate where the probability mass is, and that placement
is what buys the accuracy above at a density whose shape the integrator
knows nothing about.

## One call, not one per quantile

Every quantile is integrated in a single batched
[`numericals7::quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.html)
call, one row per quantile, so a vector of \\n\\ quantiles costs matrix
evaluations of the density in place of \\n\\ adaptive runs.

## See also

[`find_pdf_anchor()`](https://statmodels7.github.io/distributions7/reference/find_pdf_anchor.md),
which locates the mode;
[`numericals7::quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.html),
which does the integration;
[`distrib_quantile.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.continuous_distrib.md),
which inverts this;
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.
