# Skewness of the Poisson Distribution

Closed form: \\\gamma_1 = 1/\sqrt\mu\\. It is positive at every mean, a
count having a floor at zero and no ceiling, and it decays like
\\\mu^{-1/2}\\, which is the rate at which the family approaches a
Gaussian as the counts grow.

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

A numeric vector, the length of `theta$mu`, positive throughout.

## See also

[`kurtosis.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.PoissonDistrib.md),
which decays twice as fast;
[`skewness.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.NegBin2Distrib.md),
which tends to this;
[`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md).

## Examples

``` r
d <- poisson_distrib()

# One over the square root of the mean.
all.equal(skewness(d, list(mu = 4)), 0.5)
#> [1] TRUE

# It vanishes as the counts grow.
skewness(d, list(mu = c(1, 25, 10000)))
#> [1] 1.00 0.20 0.01
```
