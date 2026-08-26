# Variance of the Inverse Gaussian Distribution in Mean and Rate

Closed form: \\\operatorname{Var}(Y) = \mu^3/\lambda\\. With \\\lambda =
1/\phi\\ this is the \\\phi\mu^3\\ of
[`variance.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.InvGauss1Distrib.md),
so the two parametrizations agree on the law and differ only in which
number is reported. The cube in the mean is what separates the family
from the gamma.

## Arguments

- x:

  An `InvGauss2Distrib`, from
  [`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md).

- theta:

  A named list with components `mu` (positive) and `lambda` (positive),
  each a numeric vector of length 1 or `n`. The variance diverges as the
  rate approaches zero.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length
`max(length(theta$mu), length(theta$lambda))`.

## Notation

\\\mu \> 0\\ is the mean and \\\lambda \> 0\\ the rate.

## See also

[`mean.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.InvGauss2Distrib.md),
[`skewness.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.InvGauss2Distrib.md);
[`variance.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.InvGauss1Distrib.md),
the same quantity in the dispersion parametrization;
[`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md).

## Examples

``` r
d <- invgauss2_distrib()

# The cube of the mean over the rate.
all.equal(variance(d, list(mu = 2, lambda = 8)), 1)
#> [1] TRUE

# The two parametrizations agree when lambda = 1 / phi.
all.equal(variance(d, list(mu = 2, lambda = 8)),
          variance(invgauss1_distrib(), list(mu = 2, phi = 1 / 8)))
#> [1] TRUE
```
