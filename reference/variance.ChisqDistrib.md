# Variance of the Chi-Squared Distribution

Closed form: \\\operatorname{Var}(Y) = 2\mu\\. The variance is tied to
the mean, so the family has one parameter where a gamma has two: it is
exactly the gamma at \\\sigma^2 = 2\mu\\, and a sample whose variance is
far from twice its mean is not chi-squared at any degrees of freedom.

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

A numeric vector of variances, the length of `theta$mu`.

## See also

[`mean.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.ChisqDistrib.md);
[`variance.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Gamma2Distrib.md),
where the two moments are free of each other;
[`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md).

## Examples

``` r
d <- chisq_distrib()

# Twice the mean.
all.equal(variance(d, list(mu = 5)), 10)
#> [1] TRUE

# The gamma at sigma2 = 2 mu is the same law.
all.equal(variance(d, list(mu = 5)),
          variance(gamma2_distrib(), list(mu = 5, sigma2 = 10)))
#> [1] TRUE
```
