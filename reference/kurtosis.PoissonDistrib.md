# Excess Kurtosis of the Poisson Distribution

Closed form: \\\gamma_2 = 1/\mu\\, the excess over the Gaussian. It is
positive at every mean and decays like \\\mu^{-1}\\, twice as fast as
the skewness, so a Poisson of moderate mean looks Gaussian in its tails
before it looks symmetric.

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

A numeric vector of excess kurtoses, the length of `theta$mu`, positive
throughout.

## See also

[`skewness.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.PoissonDistrib.md);
[`kurtosis.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.NegBin2Distrib.md),
which tends to this;
[`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md).

## Examples

``` r
d <- poisson_distrib()

# One over the mean, falling twice as fast as the skewness.
rbind(skewness = skewness(d, list(mu = c(1, 25, 10000))),
      kurtosis = kurtosis(d, list(mu = c(1, 25, 10000))))
#>          [,1] [,2]  [,3]
#> skewness    1 0.20 1e-02
#> kurtosis    1 0.04 1e-04
```
