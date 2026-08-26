# Weibull Quantile Function

Computes the Weibull quantile function \$\$Q(p; \mu, \sigma) = \mu
\left\\-\log(1 - p)\right\\^{1/\sigma}\$\$ by calling
[`stats::qweibull()`](https://rdrr.io/r/stats/Weibull.html). The
distribution function is strictly increasing on \\(0, \infty)\\, so the
inverse is exact and unique; the root-finding fallback the base class
supplies is bypassed. `Q(0)` is 0 and `Q(1)` is `Inf`.

## Arguments

- distrib:

  A `Weibull1Distrib` object, from
  [`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

- p:

  A numeric vector of probabilities in \\\[0, 1\]\\, or of their
  logarithms when `log.p = TRUE`. A value outside the range gives `NaN`
  with a warning from
  [`stats::qweibull()`](https://rdrr.io/r/stats/Weibull.html).

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `p`. A component of length 1 is
  recycled. Both must be strictly positive.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, `p` is \\P(Y \le q)\\;
  when `FALSE` it is \\P(Y \> q)\\.

- log.p:

  Logical of length 1. When `TRUE`, `p` is read as a logarithm. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of quantiles in \\\[0, \infty\]\\, of length
`max(length(p), length(mu), length(sigma))`.

## See also

[`distrib_cdf.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Weibull1Distrib.md)
for the function inverted here,
[`distrib_rng.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.Weibull1Distrib.md)
for draws, and
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the generic.

## Examples

``` r
d <- weibull1_distrib()
th <- list(mu = 2, sigma = 1.5)

# The quartiles, and the round trip back through the distribution function.
q <- distrib_quantile(d, c(0.25, 0.5, 0.75), th)
q
#> [1] 0.8715759 1.5664395 2.4865678
all.equal(distrib_cdf(d, q, th), c(0.25, 0.5, 0.75))
#> [1] TRUE

# The median is mu times (log 2)^(1/sigma), whatever the shape.
all.equal(distrib_quantile(d, 0.5, th), 2 * log(2)^(1 / 1.5))
#> [1] TRUE
```
