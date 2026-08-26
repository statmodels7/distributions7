# Excess Kurtosis of the Geometric Distribution

Closed form: \\\gamma_2 = 6 + 1/\\\mu(1+\mu)\\\\, the excess over the
Gaussian. It stays above 6 at every mean and tends to 6 as the counts
grow, the exponential's value; the family is heavy-tailed at every
parameter setting and does not become Gaussian as a Poisson does.

## Arguments

- x:

  A `GeometricDistrib`, from
  [`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md).

- theta:

  A named list with one component, `mu` (positive), a numeric vector of
  length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of excess kurtoses, the length of `theta$mu`, above 6
throughout.

## See also

[`skewness.GeometricDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.GeometricDistrib.md),
which tends to 2;
[`kurtosis.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.ExponentialDistrib.md),
the continuous limit;
[`kurtosis.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.PoissonDistrib.md),
which does vanish;
[`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md).

## Examples

``` r
d <- geometric_distrib()

# The published form, written out.
all.equal(kurtosis(d, list(mu = 3)), 6 + 1 / 12)
#> [1] TRUE

# It stays above six where a Poisson's excess kurtosis vanishes.
rbind(geometric = kurtosis(d, list(mu = c(1, 10, 1000))),
      poisson   = kurtosis(poisson_distrib(), list(mu = c(1, 10, 1000))))
#>           [,1]     [,2]     [,3]
#> geometric  6.5 6.009091 6.000001
#> poisson    1.0 0.100000 0.001000
```
