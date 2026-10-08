# Student t Moments, Standard Deviation

The mean is \\\mu\\ and the variance \\\sigma^2\\, both finite on the
whole domain \\\nu \> 2\\. The skewness is zero for \\\nu \> 3\\ and
`NaN` otherwise; the excess kurtosis is \\6/(\nu - 4)\\ for \\\nu \> 4\\
and `Inf` otherwise.

## Arguments

- x:

  A `StudentT2Distrib` object.

- theta:

  A named list with components `mu`, `sigma` and `nu`.

- ...:

  Unused.

## Value

A numeric vector, one value per parameter row.

## See also

[`student_t2_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t2_distrib.md).

## Examples

``` r
d <- student_t2_distrib()
th <- list(mu = 1, sigma = 2, nu = 6)
c(mean(d, th), variance(d, th), skewness(d, th), kurtosis(d, th))
#> [1] 1 4 0 3
```
