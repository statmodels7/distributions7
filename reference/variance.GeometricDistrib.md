# Variance of the Geometric Distribution

Closed form: \\\operatorname{Var}(Y) = \mu(1+\mu)\\. It exceeds the mean
at every parameter value, so the family is overdispersed relative to a
Poisson, and it is exactly the negative binomial's variance at \\\theta
= 1\\.

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

A numeric vector of variances, the length of `theta$mu`.

## See also

[`mean.GeometricDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.GeometricDistrib.md);
[`variance.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.NegBin2Distrib.md),
of which this is the \\\theta = 1\\ case;
[`variance.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.PoissonDistrib.md),
the equidispersed comparison;
[`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md).

## Examples

``` r
d <- geometric_distrib()

# Mean times one plus the mean.
all.equal(variance(d, list(mu = 3)), 12)
#> [1] TRUE

# The negative binomial at theta = 1 is the same law.
all.equal(variance(d, list(mu = 3)),
          variance(negbin2_distrib(), list(mu = 3, theta = 1)))
#> [1] TRUE
```
