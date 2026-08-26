# Exponential Probability Density Function

Computes the exponential density \$\$f(y; \mu) =
\dfrac{1}{\mu}\exp\left(-\dfrac{y}{\mu}\right), \qquad y \ge 0,\$\$ by
calling [`stats::dexp()`](https://rdrr.io/r/stats/Exponential.html) at
`rate = 1/mu`. The parametrization here is by the mean, so the
reciprocal is taken inside the method; a reader passing `rate` where
`mu` is expected gets the reciprocal distribution.

## Arguments

- distrib:

  An `ExponentialDistrib` object, from
  [`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md).

- y:

  A numeric vector of observations. The support is \\\[0, \infty)\\;
  [`stats::dexp()`](https://rdrr.io/r/stats/Exponential.html) returns 0
  for a negative value.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `y`. A value of length 1 is recycled.
  `mu` must be strictly positive; a non-positive value gives `NaN` with
  a warning from
  [`stats::dexp()`](https://rdrr.io/r/stats/Exponential.html).

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, of length `max(length(y), length(mu))`,
one value per observation.

## See also

[`distrib_cdf.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.ExponentialDistrib.md)
for the distribution function,
[`distrib_pdf.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Gamma1Distrib.md),
which contains this at a unit shape, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- exponential_distrib()
y <- c(0.3, 1.1, 4.0)

# The method is stats::dexp at rate = 1/mu.
all.equal(distrib_pdf(d, y, list(mu = 2)), dexp(y, rate = 1 / 2))
#> [1] TRUE

# The log-density is exactly linear in y, with slope -1/mu.
diff(distrib_pdf(d, c(1, 2, 3), list(mu = 2), log = TRUE))
#> [1] -0.5 -0.5
-1 / 2
#> [1] -0.5

# The maximum is at the origin and equals 1/mu.
distrib_pdf(d, 0, list(mu = 2))
#> [1] 0.5
```
