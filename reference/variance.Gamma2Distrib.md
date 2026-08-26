# Variance of the Gamma Distribution

Closed form: \\\operatorname{Var}(Y) = \sigma^2\\. The variance is the
second parameter of this parametrization, so the method reads it off.
The shape it implies is \\a = \mu^2/\sigma^2\\, and the higher moments
are functions of it.

## Arguments

- x:

  A `Gamma2Distrib`, from
  [`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md).

- theta:

  A named list with components `mu` (the mean, positive) and `sigma2`
  (the variance, positive), each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length
`max(length(theta$mu), length(theta$sigma2))`.

## See also

[`mean.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Gamma2Distrib.md),
the other parameter;
[`skewness.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Gamma2Distrib.md)
and
[`kurtosis.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Gamma2Distrib.md),
which are functions of the implied shape;
[`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md).

## Examples

``` r
d <- gamma2_distrib()

# The second parameter is the variance.
variance(d, list(mu = 2, sigma2 = c(0.25, 1, 4)))
#> [1] 0.25 1.00 4.00

# A chi-squared with mu degrees of freedom is this family at sigma2 = 2 mu.
all.equal(variance(chisq_distrib(), list(mu = 5)),
          variance(d, list(mu = 5, sigma2 = 10)))
#> [1] TRUE
```
