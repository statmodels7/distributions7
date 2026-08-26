# Excess Kurtosis of the Exponential Distribution

Constant: \\\gamma_2 = 6\\, the excess over the Gaussian. Like the
skewness it is fixed for the whole family, there being no shape
parameter for a standardized moment to depend on.

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

A numeric vector of 6s, the length of `theta$mu`.

## See also

[`skewness.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.ExponentialDistrib.md),
the other constant;
[`kurtosis.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Gamma2Distrib.md),
which is this at shape 1;
[`kurtosis.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Weibull1Distrib.md),
which is this at shape 1 too;
[`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md).

## Examples

``` r
d <- exponential_distrib()

# Six, whatever the mean, and the two containing families agree.
c(exponential = kurtosis(d, list(mu = 3)),
  gamma       = kurtosis(gamma2_distrib(), list(mu = 3, sigma2 = 9)),
  weibull     = kurtosis(weibull1_distrib(), list(mu = 3, sigma = 1)))
#> exponential       gamma     weibull 
#>           6           6           6 
```
