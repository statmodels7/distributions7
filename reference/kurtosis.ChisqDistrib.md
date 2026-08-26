# Excess Kurtosis of the Chi-Squared Distribution

Closed form: \\\gamma_2 = 12/\mu\\, the excess over the Gaussian. It is
positive at every degrees of freedom and decays like \\\mu^{-1}\\, twice
as fast as the skewness, so a chi-squared looks Gaussian in its tails
before it looks symmetric.

## Arguments

- x:

  A `ChisqDistrib`, from
  [`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md).

- theta:

  A named list with one component, `mu` (positive), a numeric vector of
  length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of excess kurtoses, the length of `theta$mu`, positive
throughout.

## See also

[`skewness.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.ChisqDistrib.md);
[`kurtosis.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Gamma2Distrib.md)
for the containing family;
[`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md).

## Examples

``` r
d <- chisq_distrib()

# Twelve over the degrees of freedom.
all.equal(kurtosis(d, list(mu = 5)), 12 / 5)
#> [1] TRUE

# It falls twice as fast as the skewness.
rbind(skewness = skewness(d, list(mu = c(1, 10, 1000))),
      kurtosis = kurtosis(d, list(mu = c(1, 10, 1000))))
#>               [,1]      [,2]       [,3]
#> skewness  2.828427 0.8944272 0.08944272
#> kurtosis 12.000000 1.2000000 0.01200000
```
