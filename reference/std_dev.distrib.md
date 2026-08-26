# Standard Deviation of a Distribution

The square root of
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md),
evaluated at the same parameters. This is the only route to a
distribution's standard deviation: no shipped family registers a method
of its own, so whatever a family does for its variance, closed form or
quadrature, is what happens here and the root is taken after.

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
component of `theta`. `NaN` where the variance is `NaN`, `Inf` where it
is `Inf`.

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

# A Student t below two degrees of freedom has an infinite variance.
std_dev(student_t1_distrib(), list(mu = 0, sigma = 1, nu = c(1.5, 5)))
#> [1]      Inf 1.290994
```
