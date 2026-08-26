# Variance of a Distribution

Evaluates the second central moment numerically: one call of
[`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md)
at `p = 1` gives the mean, and a second at `p = 2` takes the central
moment about it. Passing the mean forward keeps the second pass to a
single quadrature. This is the fallback for a family with no closed
form, which among the shipped families means the two von Mises.

## Arguments

- x:

  A `distrib` object.

- theta:

  A named list of parameters on the parameter scale. Components may be
  vectors, giving one variance per setting.

- ...:

  Passed to
  [`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md)
  and from there to
  [`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md).

## Value

A numeric vector of variances, of length equal to the longest component
of `theta`.

## See also

[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
for the generic and the sample version,
[`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md)
for the quadrature,
[`std_dev.distrib()`](https://statmodels7.github.io/distributions7/reference/std_dev.distrib.md)
for its square root.

## Examples

``` r
# The von Mises reaches this method; on (-pi, pi] its variance is finite.
variance(vonmises1_distrib(), list(mu = 0, kappa = 2))
#> [1] 0.7644619

# Two passes of moment() give the same answer as the closed form elsewhere.
d <- gaussian1_distrib()
m <- moment(d, list(mu = 2, sigma = 3), p = 1)
all.equal(moment(d, list(mu = 2, sigma = 3), p = 2, central = TRUE, mu = m), 9)
#> [1] TRUE
```
