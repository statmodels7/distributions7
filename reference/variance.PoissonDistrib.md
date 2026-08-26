# Variance of the Poisson Distribution

Closed form: \\\operatorname{Var}(Y) = \mu\\, equal to the mean. This is
equidispersion, and it is the assumption every overdispersed count
family in the toolkit relaxes: a sample whose variance exceeds its mean
is evidence against a Poisson and for a negative binomial or a
Poisson-inverse Gaussian.

## Arguments

- x:

  A `PoissonDistrib`, from
  [`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md).

- theta:

  A named list with one component, `mu` (positive), a numeric vector of
  length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, the length of `theta$mu`.

## See also

[`mean.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.PoissonDistrib.md),
equal to this;
[`variance.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.NegBin2Distrib.md),
which adds \\\mu^2/\theta\\;
[`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md).

## Examples

``` r
d <- poisson_distrib()

# Equal to the mean at every setting.
variance(d, list(mu = c(0.5, 3, 100)))
#> [1]   0.5   3.0 100.0

# The negative binomial exceeds it and falls onto it as theta grows.
variance(negbin2_distrib(), list(mu = 3, theta = c(1, 100, 1e8)))
#> [1] 12.00  3.09  3.00
```
