# Truncated Quantile Function (Discrete)

Evaluates \\Q_T(p) = Q\\F(L^-;\theta) + p\\Z(\theta)\\\\ through the
PARENT's quantile function. The generalized inverse of a discrete
distribution function satisfies the same relation as a continuous one,
so the discrete case needs no separate treatment and returns a support
point of the truncated law.

## Arguments

- distrib:

  A `TruncatedDiscreteDistrib` object, from
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

A numeric vector of quantiles, every one a support point in \\\[L,
U\]\\.

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
[`distrib_cdf.TruncatedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.TruncatedDiscreteDistrib.md)
for the inverse, and
[`distrib_rng.TruncatedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.TruncatedDiscreteDistrib.md),
which draws through it.

## Examples

``` r
ztp <- truncated(poisson_distrib(), lower = 1)
theta <- list(mu = 2)

distrib_quantile(ztp, c(0.1, 0.5, 0.9), theta)
#> [1] 1 2 4

# Never below the lower endpoint, whatever probability is asked for.
distrib_quantile(ztp, c(0, 1e-12, 0.5), theta)
#> [1] 0 1 2
```
