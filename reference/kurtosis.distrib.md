# Excess Kurtosis of a Distribution

Evaluates the excess kurtosis numerically. Three calls of
[`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md)
are made: the mean at `p = 1`, then the second and fourth central
moments about it, and the value is \\m_4 / m_2^{2} - 3\\. Passing the
mean forward keeps the count at three quadratures. This is the fallback
for a family with no closed form, which among the shipped families means
the elastic net and the two von Mises.

## Arguments

- x:

  A `distrib` object.

- theta:

  A named list of parameters on the parameter scale. Components may be
  vectors, giving one excess kurtosis per setting.

- ...:

  Passed to
  [`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md)
  and from there to
  [`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md).

## Value

A numeric vector, of length equal to the longest component of `theta`.

## See also

[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)
for the generic and the sample version,
[`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md)
for the quadrature,
[`skewness.distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.distrib.md)
for the third order.

## Examples

``` r
# The elastic net reaches this method; its excess kurtosis lies between the
# Gaussian's 0 and the Laplace's 3.
round(kurtosis(enet_distrib(), list(mu = 0, lambda = 1, alpha = 0.5)), 4)
#> [1] 0.5523

# The same three moments, assembled by hand, on a family with a formula.
d <- gamma2_distrib(); th <- list(mu = 2, sigma2 = 1)
m <- moment(d, th, p = 1)
m2 <- moment(d, th, p = 2, central = TRUE, mu = m)
m4 <- moment(d, th, p = 4, central = TRUE, mu = m)
all.equal(m4 / m2^2 - 3, kurtosis(d, th))
#> [1] TRUE
```
