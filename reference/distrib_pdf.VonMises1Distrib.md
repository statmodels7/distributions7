# von Mises Density

Computes the von Mises density \$\$f(y; \mu, \kappa) = \dfrac{e^{\kappa
\cos(y - \mu)}}{2\pi I_0(\kappa)}, \qquad y \in \[-\pi, \pi),\$\$ with
\\I_0\\ the modified Bessel function of the first kind. Outside
\\\[-\pi, \pi)\\ the density is 0: the support is the declared interval,
and an angle is not wrapped into it.

The normalizing constant goes through
[`numericals7::log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.html).
R's exponentially scaled `besselI` **underflows to an exact zero**
between \\\kappa = 10^5\\ and \\10^6\\, where the logarithm then returns
`-Inf`; `log_bessel_i` stays finite wherever the logarithm itself is
representable. Measured at \\\kappa = 10^6\\ it returns 999992.17 where
the base route returns `-Inf`.

## Arguments

- distrib:

  A `VonMises1Distrib` object, from
  [`vonmises1_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md).

- y:

  A numeric vector of angles. A value outside \\\[-\pi, \pi)\\ is off
  the support and gives a density of 0, or `-Inf` with `log = TRUE`.

- theta:

  A named list with components `mu` and `kappa`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `mu` must lie in \\(-\pi, \pi)\\ and `kappa` be strictly
  positive.

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads
  [`numericals7::log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.html)
  may use. It is carried down because this is where the family spends
  its time: profiled at 80.8 per cent of a fit whose concentration is
  modeled, so that `kappa` is a vector. Defaults to `1L`.

## Value

A numeric vector of densities, of length
`max(length(y), length(mu), length(kappa))`, one value per observation.

## Notation

\\\mu\\ is the mean direction, \\\kappa \> 0\\ the concentration and
\\I_0\\ the modified Bessel function of the first kind of order zero,
`besselI(x, 0)` in R.

## See also

[`distrib_cdf.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.VonMises1Distrib.md)
for the distribution function,
[`distrib_gradient.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.VonMises1Distrib.md)
for the derivatives of the log-density,
[`numericals7::log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.html)
for the constant, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- vonmises1_distrib()
th <- list(mu = 0.5, kappa = 2)
y <- c(-1, 0, 0.5, 2)

# The formula written out.
all.equal(distrib_pdf(d, y, th),
          exp(2 * cos(y - 0.5)) / (2 * pi * besselI(2, 0)))
#> [1] TRUE

# It integrates to one over the circle.
integrate(function(v) distrib_pdf(d, v, th), -pi, pi)$value
#> [1] 1

# Outside the declared interval the density is zero; an angle is not
# wrapped into it.
distrib_pdf(d, c(-4, 3.5), th)
#> [1] 0 0

# The constant survives a concentration at which the base route does not.
c(ours = numericals7::log_bessel_i(1e6, 0),
  base = log(besselI(1e6, 0, expon.scaled = TRUE)) + 1e6)
#>     ours     base 
#> 999992.2     -Inf 
```
