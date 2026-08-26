# Skew Normal Density in the Centered Parametrization

Computes the skew normal density at the direct parameters that the
centered triple implies. With \\(\xi, \omega, \alpha) = \mathrm{DP}(\mu,
\sigma, \gamma_1)\\ from
[`sn_cp_to_dp()`](https://statmodels7.github.io/distributions7/reference/sn_cp_to_dp.md)
and \\z = (y-\xi)/\omega\\, \$\$f(y; \mu, \sigma, \gamma_1) =
\dfrac{2}{\omega}\\\phi(z)\\\Phi(\alpha z).\$\$ The density is the same
function of \\y\\ as
[`distrib_pdf.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.SkewNormal1Distrib.md)'s;
only the three numbers naming it differ.

Unlike the derivatives, the density is defined at \\\gamma_1 = 0\\ and
equals the Gaussian's there: the map itself is continuous through zero,
and only its derivative is not.

## Arguments

- distrib:

  A `SkewNormal2Distrib` object, from
  [`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md).

- y:

  A numeric vector of observations, anywhere on the real line.

- theta:

  A named list with components `mu`, `sigma` and `gamma1`, each a
  numeric vector of length 1 or of the length of `y`. `sigma` must be
  strictly positive and `gamma1` must lie in \\(-0.9952717,
  0.9952717)\\; a skewness outside that range belongs to no skew normal
  and the map returns `NaN`.

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, of length
`max(length(y), length(mu), length(sigma), length(gamma1))`.

## Notation

\\\mu\\, \\\sigma\\ and \\\gamma_1\\ are the mean, the standard
deviation and the skewness; \\\xi\\, \\\omega\\ and \\\alpha\\ the
location, scale and shape they imply; \\\phi\\ and \\\Phi\\ the standard
Gaussian density and distribution function.

## See also

[`distrib_pdf.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.SkewNormal1Distrib.md)
for the same density in the direct parametrization,
[`sn_cp_to_dp()`](https://statmodels7.github.io/distributions7/reference/sn_cp_to_dp.md)
for the map, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- skewnormal2_distrib()
d1 <- skewnormal1_distrib()
y <- c(-1, 0.3, 1.7)
th <- list(mu = 0, sigma = 1, gamma1 = 0.5)

# The same law, reached through the map.
all.equal(distrib_pdf(d, y, th),
          distrib_pdf(d1, y, distributions7:::sn2_theta(th)))
#> [1] TRUE

# It integrates to one, and its first moment is the parameter mu.
c(mass = integrate(function(v) distrib_pdf(d, v, th), -Inf, Inf)$value,
  mean = integrate(function(v) v * distrib_pdf(d, v, th), -Inf, Inf)$value)
#>          mass          mean 
#>  1.000000e+00 -2.257916e-14 

# At zero skewness the density is the Gaussian's, where the derivatives
# are not defined.
all.equal(distrib_pdf(d, y, list(mu = 0, sigma = 1, gamma1 = 0)), dnorm(y))
#> [1] TRUE
```
