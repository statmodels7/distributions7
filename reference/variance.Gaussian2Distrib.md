# Variance of the Gaussian Distribution in Location and Variance

Closed form: \\\operatorname{Var}(Y) = \sigma^2\\. The variance is the
second parameter here, so the method reads it off, where
[`variance.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Gaussian1Distrib.md)
squares a standard deviation and
[`variance.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Gaussian3Distrib.md)
inverts a precision. All three describe the same law.

## Arguments

- x:

  A `Gaussian2Distrib`, from
  [`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md).

- theta:

  A named list with components `mu` (any real value) and `sigma2` (the
  variance, positive), each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length
`max(length(theta$mu), length(theta$sigma2))`.

## See also

[`mean.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Gaussian2Distrib.md);
[`variance.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Gaussian1Distrib.md)
and
[`variance.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Gaussian3Distrib.md),
the same quantity in the other two parametrizations;
[`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md).

## Examples

``` r
d <- gaussian2_distrib()

# The second parameter is the variance.
variance(d, list(mu = 0, sigma2 = c(1, 4, 16)))
#> [1]  1  4 16

# The three parametrizations agree on the law.
c(gaussian2 = variance(d, list(mu = 0, sigma2 = 9)),
  gaussian1 = variance(gaussian1_distrib(), list(mu = 0, sigma = 3)),
  gaussian3 = variance(gaussian3_distrib(), list(mu = 0, tau = 1 / 9)))
#> gaussian2 gaussian1 gaussian3 
#>         9         9         9 
```
