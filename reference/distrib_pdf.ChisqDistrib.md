# Chi-Squared Probability Density Function

Computes the chi-squared density \$\$f(y; \mu) = \dfrac{y^{\mu/2 - 1}
e^{-y/2}}{2^{\mu/2}\\\Gamma(\mu/2)}, \qquad y \> 0,\$\$ by calling
[`stats::dchisq()`](https://rdrr.io/r/stats/Chisquare.html) at
`df = mu`. With `log = TRUE` the logarithm is formed inside
[`dchisq()`](https://rdrr.io/r/stats/Chisquare.html) and stays finite
where the density itself underflows.

The density is unbounded at the origin for \\\mu \< 2\\, flat there at
\\\mu = 2\\, where the family is the exponential with mean 2, and
vanishes at the origin for \\\mu \> 2\\.

## Arguments

- distrib:

  A `ChisqDistrib` object, from
  [`chisq_distrib()`](https://statmodels7.github.io/distributions7/reference/chisq_distrib.md).

- y:

  A numeric vector of observations. The support is \\(0, \infty)\\; a
  value at or below zero gives 0, `Inf` or a finite value according to
  \\\mu\\.

- theta:

  A named list with one component `mu`, a numeric vector of length 1 or
  of the length of `y`, recycled if of length 1. It must be strictly
  positive and need not be a whole number.

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, of length `max(length(y), length(mu))`,
one value per observation.

## See also

[`distrib_cdf.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.ChisqDistrib.md)
for the distribution function,
[`distrib_gradient.ChisqDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ChisqDistrib.md)
for the derivative of the log-density,
[`distrib_pdf.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Gamma2Distrib.md)
for the family this sits inside, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- chisq_distrib()
y <- c(1, 4, 9)
th <- list(mu = 4)

# The method is stats::dchisq at df = mu.
all.equal(distrib_pdf(d, y, th), dchisq(y, df = 4))
#> [1] TRUE

# It is the gamma with sigma2 = 2 mu, and the exponential at mu = 2.
all.equal(distrib_pdf(d, y, th),
          distrib_pdf(gamma2_distrib(), y, list(mu = 4, sigma2 = 8)))
#> [1] TRUE
all.equal(distrib_pdf(d, y, list(mu = 2)),
          distrib_pdf(exponential_distrib(), y, list(mu = 2)))
#> [1] TRUE

# The degrees of freedom need not be a whole number.
distrib_pdf(d, y, list(mu = 3.7))
#> [1] 0.17792406 0.12898641 0.02109419

# Far out in the tail the density underflows and its logarithm does not.
distrib_pdf(d, 2000, th)
#> [1] 0
distrib_pdf(d, 2000, th, log = TRUE)
#> [1] -993.7854
```
