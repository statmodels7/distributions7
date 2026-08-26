# Student t Probability Density Function

Computes the location-scale Student t density \$\$f(y; \mu, \sigma, \nu)
=
\dfrac{\Gamma\left(\dfrac{\nu+1}{2}\right)}{\sigma\sqrt{\nu\pi}\\\Gamma\left(\dfrac{\nu}{2}\right)}
\left(1 + \dfrac{(y-\mu)^2}{\nu\sigma^2}\right)^{-\dfrac{\nu+1}{2}}\$\$
by calling [`stats::dt()`](https://rdrr.io/r/stats/TDist.html) at the
standardized value \\(y-\mu)/\sigma\\ and dividing by \\\sigma\\, the
Jacobian of the standardization. With `log = TRUE` the division becomes
a subtraction of \\\log\sigma\\, so the logarithm stays finite far into
the tails.

The tails are polynomial, of order \\\|y\|^{-(\nu+1)}\\, so the density
decays far more slowly than a Gaussian's and no value of \\y\\
underflows at an ordinary \\\nu\\.

## Arguments

- distrib:

  A `StudentT1Distrib` object, from
  [`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

- y:

  A numeric vector of observations. Every real value is in the support,
  so no value is rejected.

- theta:

  A named list with components `mu`, `sigma` and `nu`, each a numeric
  vector of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` and `nu` must be strictly positive; a zero or
  negative value gives `NaN` with a warning from
  [`stats::dt()`](https://rdrr.io/r/stats/TDist.html).

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, of length
`max(length(y), length(mu), length(sigma), length(nu))`, one value per
observation.

## See also

[`distrib_cdf.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.StudentT1Distrib.md)
for the distribution function,
[`distrib_gradient.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.StudentT1Distrib.md)
for the derivatives of the log-density,
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic and
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md)
for the family.

## Examples

``` r
d <- student_t1_distrib()
y <- c(-2.5, 0.3, 1.8)

# The method is stats::dt at the standardized value, over sigma.
all.equal(distrib_pdf(d, y, list(mu = 0.4, sigma = 1.2, nu = 5)),
          dt((y - 0.4) / 1.2, df = 5) / 1.2)
#> [1] TRUE

# One degree of freedom is the Cauchy.
all.equal(distrib_pdf(d, y, list(mu = 0.4, sigma = 1.2, nu = 1)),
          dcauchy(y, location = 0.4, scale = 1.2))
#> [1] TRUE

# The tails are polynomial, so a value forty scales out still carries mass
# where a Gaussian's density has underflowed to zero.
c(t = distrib_pdf(d, 48, list(mu = 0, sigma = 1.2, nu = 5)),
  gaussian = dnorm(48, 0, 1.2))
#>            t     gaussian 
#> 9.563955e-09 0.000000e+00 
```
