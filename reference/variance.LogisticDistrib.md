# Variance of the Logistic Distribution

Closed form: \\\operatorname{Var}(Y) = \pi^2\sigma^2/3\\. The standard
deviation is \\\sigma\pi/\sqrt3 \approx 1.8138\sigma\\, so the scale is
smaller than the standard deviation by the factor \\\pi/\sqrt3\\, and a
logistic scale is not comparable with a Gaussian standard deviation as
it stands.

## Arguments

- x:

  A `LogisticDistrib`, from
  [`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md).

- theta:

  A named list with components `mu` (the location) and `sigma` (the
  scale, positive), each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length
`max(length(theta$mu), length(theta$sigma))`.

## See also

[`mean.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.LogisticDistrib.md),
[`kurtosis.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.LogisticDistrib.md),
[`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md).

## Examples

``` r
d <- logistic_distrib()

# pi^2 sigma^2 / 3.
all.equal(variance(d, list(mu = 0, sigma = 2)), pi^2 * 4 / 3)
#> [1] TRUE

# The standard deviation is pi / sqrt(3) times the scale.
std_dev(d, list(mu = 0, sigma = 1))
#> [1] 1.813799
```
