# Truncated Cumulative Distribution Function (Continuous)

Evaluates \\F_T(q) = \\F(q;\theta) - F(L^-;\theta)\\/Z(\theta)\\,
clamped to \\\[0, 1\]\\. The clamp makes the endpoints exact, returning
`0` at and below \\L\\ and `1` at and above \\U\\, where the unclamped
ratio would give a small negative number or one plus a rounding.

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object, from
  [`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md).

- q:

  A numeric vector of quantiles.

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
[`distrib_quantile.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.TruncatedContinuousDistrib.md)
for the inverse, and
[`distrib_cdf.TruncatedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.TruncatedDiscreteDistrib.md)
for the discrete branch.

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

# Exactly 0 and 1 at the endpoints.
distrib_cdf(tn, c(-1, 0.3, 2), theta)
#> [1] 0.0000000 0.4609908 1.0000000

# Against the ratio written out.
Fl <- pnorm(-1, 0.3, 1.2); Z <- pnorm(2, 0.3, 1.2) - Fl
(pnorm(0.3, 0.3, 1.2) - Fl) / Z
#> [1] 0.4609908

# The two tails sum to one.
distrib_cdf(tn, 0.3, theta) + distrib_cdf(tn, 0.3, theta, lower.tail = FALSE)
#> [1] 1
```
