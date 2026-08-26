# Variance of the Skew Normal in the Centered Parametrization

Returns \\\sigma^2\\, the square of the second parameter, which in this
parametrization is the variance by construction.
[`moment_const()`](https://statmodels7.github.io/distributions7/reference/moment_const.md)
recycles the result to the length the three parameters imply.

## Arguments

- x:

  A `SkewNormal2Distrib` object, from
  [`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md).

- theta:

  A named list with components `mu`, `sigma` and `gamma1`, in any order;
  it is aligned here.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of the length the recycled parameters
imply.

## See also

[`mean.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.SkewNormal2Distrib.md)
and
[`skewness.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.SkewNormal2Distrib.md),
and
[`variance.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.SkewNormal1Distrib.md),
where the scale is not the standard deviation.

## Examples

``` r
d <- skewnormal2_distrib()
variance(d, list(mu = 3, sigma = 2, gamma1 = 0.6))
#> [1] 4

# It agrees with a quadrature of the second central moment.
th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
c(parameter = variance(d, th),
  quadrature = integrate(function(v) v^2 * distrib_pdf(d, v, th),
                         -Inf, Inf)$value)
#>  parameter quadrature 
#>          1          1 
```
