# Standard Deviation of a Distribution

Returns the square root of
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md),
evaluated at the same parameters. This is the only route to a
distribution's standard deviation, because a shipped family registers no
method of its own. The variance of the family, by closed form or by
quadrature, is therefore computed first and the root is taken after.

## Arguments

- x:

  A `distrib` object.

- theta:

  A named list of parameters on the parameter scale. Components may be
  vectors, giving one standard deviation per setting.

- ...:

  Passed to
  [`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
  and from there to
  [`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md).

## Value

A numeric vector of standard deviations, of length equal to the longest
component of `theta`. The value is `NaN` where the variance is `NaN` and
`Inf` where it is `Inf`.

## See also

[`std_dev()`](https://statmodels7.github.io/distributions7/reference/std_dev.md)
for the generic and the sample version,
[`variance.distrib()`](https://statmodels7.github.io/distributions7/reference/variance.distrib.md)
for the quantity under the root.

## Examples

``` r
d <- weibull1_distrib()
th <- list(mu = 2, sigma = 3)
all.equal(std_dev(d, th), sqrt(variance(d, th)))
#> [1] TRUE

# A Student t has an infinite variance for 1 < nu <= 2 (NaN for nu <= 1).
std_dev(student_t1_distrib(), list(mu = 0, sigma = 1, nu = c(1.5, 5)))
#> [1]      Inf 1.290994
```
