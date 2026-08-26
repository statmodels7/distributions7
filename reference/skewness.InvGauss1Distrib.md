# Skewness of the Inverse Gaussian Distribution

Closed form: \\\gamma_1 = 3\sqrt{\phi\mu}\\. It is positive at every
parameter value and depends on the two parameters only through their
product, so the whole family is one curve indexed by \\\phi\mu\\.

## Arguments

- x:

  An `InvGauss1Distrib`, from
  [`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md).

- theta:

  A named list with components `mu` (positive) and `phi` (positive),
  each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector, of length `max(length(theta$mu), length(theta$phi))`,
positive throughout.

## Notation

\\\mu \> 0\\ is the mean and \\\phi \> 0\\ the dispersion.

## See also

[`kurtosis.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.InvGauss1Distrib.md),
which is \\5/3\\ times the square of this;
[`variance.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.InvGauss1Distrib.md);
[`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md).

## Examples

``` r
d <- invgauss1_distrib()

# Three times the square root of the product.
all.equal(skewness(d, list(mu = 2, phi = 0.5)), 3 * sqrt(1))
#> [1] TRUE

# Only the product matters.
skewness(d, list(mu = c(1, 2, 4), phi = c(4, 2, 1)))
#> [1] 6 6 6
```
