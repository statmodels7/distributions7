# Poisson Probability Mass Function

Computes the Poisson probability mass \$\$P(Y = y; \mu) =
\dfrac{e^{-\mu}\mu^{y}}{y!}, \qquad y = 0, 1, 2, \dots\$\$ by calling
[`stats::dpois()`](https://rdrr.io/r/stats/Poisson.html). The `pdf` in
the generic's name is the density with respect to counting measure, so
what this returns is a probability and is at most 1.

## Arguments

- distrib:

  A `PoissonDistrib` object, from
  [`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md).

- y:

  A numeric vector of counts. A non-integer or negative value gives 0
  with a warning from
  [`stats::dpois()`](https://rdrr.io/r/stats/Poisson.html).

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `y`. A value of length 1 is recycled.
  `mu` must be strictly positive; a negative value gives `NaN` with a
  warning.

- log:

  Logical of length 1. When `TRUE` the log-mass is returned, which stays
  finite for a count far above the mean where the mass itself
  underflows. Defaults to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, of length
`max(length(y), length(mu))`, one value per observation.

## See also

[`distrib_cdf.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.PoissonDistrib.md)
for the distribution function,
[`distrib_gradient.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.PoissonDistrib.md)
for the derivatives of the log-mass, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- poisson_distrib()
y <- c(0, 2, 7)

# The method is stats::dpois.
all.equal(distrib_pdf(d, y, list(mu = 3)), dpois(y, lambda = 3))
#> [1] TRUE

# A probability mass: it sums to one over the support.
sum(distrib_pdf(d, 0:200, list(mu = 3)))
#> [1] 1

# Far above the mean the mass underflows and its logarithm does not.
distrib_pdf(d, 400, list(mu = 3))
#> [1] 0
distrib_pdf(d, 400, list(mu = 3), log = TRUE)
#> [1] -1564.056
```
