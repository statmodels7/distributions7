# Variance of the Laplace Distribution in Location and Rate

Closed form, replacing the numerical default: \\\operatorname{Var}(Y) =
2/\lambda^2\\. With \\\lambda = 1/\sigma\\ this is the \\2\sigma^2\\ of
[`variance.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.LaplaceDistrib.md),
so the two parametrizations agree on the law and differ only in which
number is reported. A larger rate is a tighter distribution.

## Arguments

- x:

  A `Laplace2Distrib`, from
  [`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md).

- theta:

  A named list with components `mu` (the location) and `lambda` (the
  rate, positive), each a numeric vector of length 1 or `n`. The
  variance diverges as the rate approaches zero.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length
`max(length(theta$mu), length(theta$lambda))`. The location does not
enter the value, so a setting that varies `mu` alone repeats one number.

## See also

[`variance.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.LaplaceDistrib.md),
the same quantity in the scale parametrization;
[`mean.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Laplace2Distrib.md);
[`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md).

## Examples

``` r
d <- laplace2_distrib()

# Two over the square of the rate.
all.equal(variance(d, list(mu = 0, lambda = 0.5)), 8)
#> [1] TRUE

# The two parametrizations agree when lambda = 1 / sigma.
all.equal(variance(d, list(mu = 0, lambda = 1 / 3)),
          variance(laplace_distrib(), list(mu = 0, sigma = 3)))
#> [1] TRUE
```
