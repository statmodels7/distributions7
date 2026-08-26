# Skewness of a Distribution

Evaluates the third standardized central moment numerically. Three calls
of
[`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md)
are made: the mean at `p = 1`, then the second and third central moments
about it, and the value is \\m_3 / m_2^{3/2}\\. Passing the mean forward
keeps the count at three quadratures. This is the fallback for a family
with no closed form, which among the shipped families means the two von
Mises.

## Arguments

- x:

  A `distrib` object.

- theta:

  A named list of parameters on the parameter scale. Components may be
  vectors, giving one skewness per setting.

- ...:

  Passed to
  [`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md)
  and from there to
  [`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md).

## Value

A numeric vector, of length equal to the longest component of `theta`.

## See also

[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
for the generic and the sample version,
[`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md)
for the quadrature,
[`kurtosis.distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.distrib.md)
for the fourth order.

## Examples

``` r
# The von Mises reaches this method, and is symmetric about mu.
round(skewness(vonmises1_distrib(), list(mu = 0, kappa = 2)), 10)
#> [1] 0

# The same three moments, assembled by hand, on a family that has a formula.
d <- gamma2_distrib(); th <- list(mu = 2, sigma2 = 1)
m <- moment(d, th, p = 1)
m2 <- moment(d, th, p = 2, central = TRUE, mu = m)
m3 <- moment(d, th, p = 3, central = TRUE, mu = m)
all.equal(m3 / m2^1.5, skewness(d, th))
#> [1] TRUE
```
