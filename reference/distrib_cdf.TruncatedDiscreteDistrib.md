# Truncated Cumulative Distribution Function (Discrete)

Evaluates \\F_T(q) = \\F(q;\theta) - F(L^-;\theta)\\/Z(\theta)\\,
clamped to \\\[0, 1\]\\. The lower tail subtracted is \\F(L) - f(L)\\,
the mass at \\L\\ being retained, so the truncated distribution function
is already positive at \\L\\ itself.

## Arguments

- distrib:

  A `TruncatedDiscreteDistrib` object, from
  [`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md).

- q:

  A numeric vector of quantiles, not necessarily support points.

- theta:

  A named list of the parent's parameters.

- lower.tail:

  Logical, default `TRUE`. When `FALSE` the upper tail \\1 - F_T(q)\\ is
  returned.

- log.p:

  Logical, default `FALSE`. When `TRUE` the probability is returned on
  the log scale.

## Value

A numeric vector of cumulative probabilities.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md),
[`trunc_cdf()`](https://statmodels7.github.io/distributions7/reference/trunc_cdf.md)
for the shared body,
[`distrib_quantile.TruncatedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.TruncatedDiscreteDistrib.md)
for the inverse, and
[`distrib_cdf.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.TruncatedContinuousDistrib.md)
for the continuous branch.

## Examples

``` r
ztp <- truncated(poisson_distrib(), lower = 1)
theta <- list(mu = 2)

# Zero below the interval, and already positive at the lower endpoint.
distrib_cdf(ztp, 0:4, theta)
#> [1] 0.0000000 0.3130353 0.6260706 0.8347608 0.9391059

# The cumulative sums of the truncated mass function.
cumsum(distrib_pdf(ztp, 1:4, theta))
#> [1] 0.3130353 0.6260706 0.8347608 0.9391059
```
