# von Mises Density in the Resultant Length

Computes the von Mises density \$\$f(y; \mu, \rho) =
\dfrac{e^{\kappa\cos(y-\mu)}}{2\pi I_0(\kappa)}, \qquad \kappa =
A^{-1}(\rho),\$\$ by inverting the map once and calling
[`distrib_pdf.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.VonMises1Distrib.md)
at the implied concentration. The two parametrizations are the same law,
so the value is identical to that family's at \\\kappa\\.

## Arguments

- distrib:

  A `VonMises2Distrib` object, from
  [`vonmises2_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises2_distrib.md).

- y:

  A numeric vector of angles. A value outside \\\[-\pi, \pi)\\ is off
  the support and gives a density of 0.

- theta:

  A named list with components `mu` and `rho`, each a numeric vector of
  length 1 or of the length of `y`. `mu` must lie in \\(-\pi, \pi)\\ and
  `rho` in \\(0, 1)\\.

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, passed on to
  [`numericals7::log_bessel_i()`](https://statmodels7.github.io/numericals7/reference/log_bessel_i.html)
  through the concentration parametrization. Defaults to `1L`.

## Value

A numeric vector of densities, of length
`max(length(y), length(mu), length(rho))`, one value per observation.

## Notation

\\\mu\\ is the mean direction, \\\rho \in (0,1)\\ the mean resultant
length, \\\kappa\\ the concentration, \\I_m\\ the modified Bessel
function of the first kind of order \\m\\, and \\A(\kappa) =
I_1(\kappa)/I_0(\kappa)\\.

## See also

[`distrib_pdf.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.VonMises1Distrib.md),
which this calls;
[`distrib_cdf.VonMises2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.VonMises2Distrib.md)
for the distribution function; and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d2 <- vonmises2_distrib()
th <- list(mu = 0.5, rho = 0.7)
y <- c(-1, 0, 0.5, 2)
distrib_pdf(d2, y, th)
#> [1] 0.07974218 0.40483719 0.51800669 0.07974218

# It is the same law as the concentration parametrization, at the
# concentration this resultant length implies.
k <- numericals7::bessel_i_ratio_inverse(0.7)$kappa
k
#> [1] 2.013628
all.equal(distrib_pdf(d2, y, th),
          distrib_pdf(vonmises1_distrib(), y, list(mu = 0.5, kappa = k)))
#> [1] TRUE

# And that concentration maps back to the resultant length given.
numericals7::bessel_i_ratio(k)
#> [1] 0.7
```
