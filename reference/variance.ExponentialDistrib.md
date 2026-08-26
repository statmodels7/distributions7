# Variance of the Exponential Distribution

Closed form: \\\operatorname{Var}(Y) = \mu^2\\. The standard deviation
equals the mean, so the coefficient of variation is exactly one at every
parameter value. That single number is what identifies an exponential
among the gammas: a sample whose relative spread is far from one is not
one.

## Arguments

- x:

  An `ExponentialDistrib`, from
  [`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md).

- theta:

  A named list with one component, `mu` (the mean, positive), a numeric
  vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, the length of `theta$mu`.

## See also

[`mean.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.ExponentialDistrib.md);
[`variance.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Gamma2Distrib.md),
where the two moments are free of each other;
[`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md).

## Examples

``` r
d <- exponential_distrib()

# The square of the mean, so the coefficient of variation is one.
variance(d, list(mu = 3))
#> [1] 9
std_dev(d, list(mu = 3)) / mean(d, list(mu = 3))
#> [1] 1

# The same law as a Weibull of shape 1.
all.equal(variance(d, list(mu = 3)),
          variance(weibull1_distrib(), list(mu = 3, sigma = 1)))
#> [1] TRUE
```
