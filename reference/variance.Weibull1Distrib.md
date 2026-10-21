# Variance of the Weibull Distribution

The closed form is \\\operatorname{Var}(Y) = \mu^2 (g_2 - g_1^2)\\ with
\\g_k = \Gamma(1 + k/\sigma)\\. The scale enters as a square and the
shape through the two gamma factors, so the coefficient of variation
\\\sqrt{g_2 - g_1^2}/g_1\\ depends only on the shape. The method
computes \\\mu^2 g_1^2 (g_2/g_1^2 - 1)\\ from the logarithms of the
gamma factors
([`weibull_moment_terms()`](https://statmodels7.github.io/distributions7/reference/weibull_moment_terms.md)),
so that \\g_2 - g_1^2\\ is not formed as a difference; the value is
`Inf` only where it exceeds the largest double, as at a shape of 0.01.

## Arguments

- x:

  A `Weibull1Distrib`, from
  [`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

- theta:

  A named list with components `mu` (the scale, positive) and `sigma`
  (the shape, positive); each is a numeric vector of length 1 or `n`.
  The variance exists for every positive shape and grows without bound
  as the shape tends to zero.

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
[`weibull_moment_terms()`](https://statmodels7.github.io/distributions7/reference/weibull_moment_terms.md),
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

## Examples

``` r
d <- weibull1_distrib()

# The variance is written out with the two gamma factors.
g1 <- gamma(1 + 1 / 3); g2 <- gamma(1 + 2 / 3)
all.equal(variance(d, list(mu = 2, sigma = 3)), 4 * (g2 - g1^2))
#> [1] TRUE

# At shape 1 the variance is the square of the scale, as for an exponential.
variance(d, list(mu = 2, sigma = 1))
#> [1] 4
```
