# Skewness of the Chi-Squared Distribution

Closed form: \\\gamma_1 = \sqrt{8/\mu}\\. It is positive at every
degrees of freedom, the support being the positive half-line, and it
decays like \\\mu^{-1/2}\\, which is the rate at which the family
approaches a Gaussian as its degrees of freedom grow.

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

A numeric vector, the length of `theta$mu`, positive throughout.

## See also

[`kurtosis.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.ChisqDistrib.md),
which decays twice as fast;
[`skewness.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Gamma2Distrib.md)
for the containing family;
[`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md).

## Examples

``` r
d <- chisq_distrib()

# Square root of eight over the degrees of freedom.
all.equal(skewness(d, list(mu = 5)), sqrt(8 / 5))
#> [1] TRUE

# It vanishes as the degrees of freedom grow.
skewness(d, list(mu = c(1, 10, 1000)))
#> [1] 2.82842712 0.89442719 0.08944272
```
