# Variance of the Weibull Distribution

Closed form: \\\operatorname{Var}(Y) = \mu^2 (g_2 - g_1^2)\\ with \\g_k
= \Gamma(1 + k/\sigma)\\. The scale enters as a square and the shape
through the two gamma factors, so the coefficient of variation
\\\sqrt{g_2 - g_1^2}/g_1\\ is a function of the shape alone.

## Arguments

- x:

  A `Weibull1Distrib`, from
  [`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

- theta:

  A named list with components `mu` (the scale, positive) and `sigma`
  (the shape, positive), each a numeric vector of length 1 or `n`. The
  variance diverges as the shape falls below \\1/2\\, the second moment
  failing to exist there.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length
`max(length(theta$mu), length(theta$sigma))`.

## Notation

\\\mu \> 0\\ is the scale, \\\sigma \> 0\\ the shape and \\g_k =
\Gamma(1 + k/\sigma)\\.

## See also

[`mean.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Weibull1Distrib.md),
[`skewness.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Weibull1Distrib.md),
[`weibull_gamma_factors()`](https://statmodels7.github.io/distributions7/reference/weibull_gamma_factors.md),
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

## Examples

``` r
d <- weibull1_distrib()

# The two gamma factors, written out.
g1 <- gamma(1 + 1 / 3); g2 <- gamma(1 + 2 / 3)
all.equal(variance(d, list(mu = 2, sigma = 3)), 4 * (g2 - g1^2))
#> [1] TRUE

# At shape 1 the variance is the square of the scale, as for an exponential.
variance(d, list(mu = 2, sigma = 1))
#> [1] 4
```
