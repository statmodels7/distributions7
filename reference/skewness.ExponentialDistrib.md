# Skewness of the Exponential Distribution

Constant: \\\gamma_1 = 2\\. The family is a scale family with no shape
parameter, so a standardized moment is a number: every exponential has
the same asymmetry, whatever its mean.

## Arguments

- x:

  An `ExponentialDistrib`, from
  [`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md).

- theta:

  A named list with one component, `mu` (positive), a numeric vector of
  length 1 or `n`. The value is not read, only its length.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of 2s, the length of `theta$mu`.

## See also

[`kurtosis.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.ExponentialDistrib.md),
the other constant;
[`skewness.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Gamma2Distrib.md),
which is this at shape 1;
[`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md).

## Examples

``` r
d <- exponential_distrib()

# Two, whatever the mean.
skewness(d, list(mu = c(0.1, 3, 100)))
#> [1] 2 2 2
```
