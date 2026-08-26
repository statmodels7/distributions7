# Variance of the Gaussian Distribution in Location and Precision

Closed form: \\\operatorname{Var}(Y) = 1/\tau\\. Larger precision means
a tighter distribution, and the name records that. The reciprocal is the
whole of the difference between this parametrization and
[`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md).

## Arguments

- x:

  A `Gaussian3Distrib`, from
  [`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md).

- theta:

  A named list with components `mu` (any real value) and `tau` (the
  precision, positive), each a numeric vector of length 1 or `n`. The
  variance diverges as the precision approaches zero.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length
`max(length(theta$mu), length(theta$tau))`.

## See also

[`mean.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Gaussian3Distrib.md);
[`variance.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Gaussian2Distrib.md),
the reciprocal;
[`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md).

## Examples

``` r
d <- gaussian3_distrib()

# One over the precision.
all.equal(variance(d, list(mu = 0, tau = 0.25)), 4)
#> [1] TRUE

# The standard deviation is one over its square root.
std_dev(d, list(mu = 0, tau = 4))
#> [1] 0.5
```
