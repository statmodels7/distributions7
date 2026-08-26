# Transformed Quantile Function

Computes \\Q_Y(p) = g(Q_X(p))\\ by evaluating the parent's own quantile
function and mapping the result forward. For a DECREASING transformation
the tails swap, so `lower.tail` is inverted before the parent is called
and the \\p\\-th quantile of \\Y\\ is \\g\\ of the \\(1-p)\\-th of
\\X\\.

## Arguments

- distrib:

  A `TransformedDistrib` object, from
  [`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md).

- p:

  A numeric vector of probabilities.

- theta:

  A named list of the parent's parameters.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\.
  For a decreasing transformation this is inverted before reaching the
  parent.

- log.p:

  Logical of length 1. When `TRUE`, `p` is given as a logarithm and
  passed as such to the parent. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of quantiles on the transformed scale.

## Notation

\\g\\ is the transformation and \\Q_X\\, \\Q_Y\\ the parent's and the
transformed quantile function.

## See also

[`distrib_cdf.TransformedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.TransformedDistrib.md),
which this inverts, and
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- transformation(gaussian1_distrib(), exp_transform())
theta <- list(mu = 0.5, sigma = 0.8)
p <- c(0.1, 0.5, 0.9)

distrib_quantile(d, p, theta)
#> [1] 0.5914127 1.6487213 4.5962523
all.equal(distrib_quantile(d, p, theta), qlnorm(p, 0.5, 0.8))
#> [1] TRUE

# The round trip closes, the transformation being a bijection.
max(abs(distrib_cdf(d, distrib_quantile(d, p, theta), theta) - p))
#> [1] 5.551115e-17

# Under a decreasing transformation the p-th quantile is g of the
# (1 - p)-th of the parent.
ig <- transformation(gamma1_distrib(), inverse_transform())
th2 <- list(mu = 2, phi = 0.3)
c(transformed = distrib_quantile(ig, 0.9, th2),
  mapped = 1 / distrib_quantile(gamma1_distrib(), 0.1, th2))
#> transformed      mapped 
#>    1.272226    1.272226 
```
