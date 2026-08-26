# Variance of the Beta Distribution

Closed form: \\\operatorname{Var}(Y) = \mu(1-\mu)/(\phi+1)\\. The
numerator is the Bernoulli variance at the same mean, the largest a
distribution on \\(0,1)\\ can have, and \\\phi+1\\ divides it down.
Larger \\\phi\\ therefore means a tighter distribution, which is the
reason the parameter is called a precision.

## Arguments

- x:

  A `Beta1Distrib`, from
  [`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md).

- theta:

  A named list with components `mu` (strictly between 0 and 1) and `phi`
  (positive), each a numeric vector of length 1 or `n`. The variance
  tends to \\\mu(1-\mu)\\ as the precision goes to zero, where the mass
  piles up at the two endpoints.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length
`max(length(theta$mu), length(theta$phi))`.

## Notation

\\\mu \in (0,1)\\ is the mean and \\\phi \> 0\\ the precision.

## See also

[`mean.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Beta1Distrib.md),
[`skewness.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Beta1Distrib.md),
[`kurtosis.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Beta1Distrib.md),
[`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md).

## Examples

``` r
d <- beta1_distrib()

# mu (1 - mu) / (phi + 1).
all.equal(variance(d, list(mu = 0.3, phi = 5)), 0.3 * 0.7 / 6)
#> [1] TRUE

# At mu = 1/2 and phi = 2 the density is uniform, of variance 1/12.
all.equal(variance(d, list(mu = 0.5, phi = 2)), 1 / 12)
#> [1] TRUE
```
