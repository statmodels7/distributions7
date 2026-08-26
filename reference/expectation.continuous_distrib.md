# Expectation of a Continuous Distribution

Evaluates \\E\[f(Y)\] = \int f(y)\\p(y;\theta)\\dy\\ by the batched
adaptive quadrature of
[`numericals7::quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.html):
the panels of every parameter combination are refined in one call, so a
vector `theta` costs matrix evaluations rather than one adaptive run per
value. The domain of each combination is split at its 0.1, 0.5 and 0.9
quantiles, which anchors the quadrature on the probability mass wherever
it sits. A combination the batched quadrature rejects – an integrable
endpoint singularity too harsh for bisection – is rescued by one scalar
[`stats::integrate()`](https://rdrr.io/r/stats/integrate.html) run,
whose extrapolation reaches it; an error naming the combination is
raised only when both routes fail.

## Arguments

- distrib:

  A `continuous_distrib`.

- f:

  The function whose expectation is taken.

- theta:

  A named list of parameters.

- ...:

  Further arguments passed to `f`.

## Value

A numeric vector of expected values.
