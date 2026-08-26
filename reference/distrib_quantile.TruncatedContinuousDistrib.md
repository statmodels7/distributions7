# Truncated Quantile Function (Continuous)

Evaluates \\Q_T(p) = Q\\F(L^-;\theta) + p\\Z(\theta)\\\\, inverting
[`distrib_cdf.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.TruncatedContinuousDistrib.md)
through the PARENT's quantile function. No root-finding of its own is
needed, and every returned value lies inside \\\[L, U\]\\ by
construction.

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object, from
  [`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md).

- p:

  A numeric vector of probabilities, clamped to \\\[0, 1\]\\.

- theta:

  A named list of the parent's parameters.

- lower.tail:

  Logical, default `TRUE`. When `FALSE`, `p` is read as an upper-tail
  probability.

- log.p:

  Logical, default `FALSE`. When `TRUE`, `p` is given on the log scale.

## Value

A numeric vector of quantiles.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md),
[`trunc_quantile()`](https://statmodels7.github.io/distributions7/reference/trunc_quantile.md)
for the shared body,
[`distrib_cdf.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.TruncatedContinuousDistrib.md)
for the inverse, and
[`distrib_rng.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.TruncatedContinuousDistrib.md),
which draws through it.

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

distrib_quantile(tn, c(0.1, 0.5, 0.9), theta)
#> [1] -0.6365189  0.3918926  1.5105956

# It inverts the distribution function exactly.
p <- c(0.1, 0.5, 0.9)
max(abs(distrib_cdf(tn, distrib_quantile(tn, p, theta), theta) - p))
#> [1] 2.220446e-16

# The extreme quantiles are the endpoints themselves.
distrib_quantile(tn, c(0, 1), theta)
#> [1] -1  2
```
