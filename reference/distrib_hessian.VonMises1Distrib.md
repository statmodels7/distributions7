# von Mises Observed Hessian

Computes the three distinct second derivatives of the von Mises
log-density, one value per observation, in closed form:
\$\$\ell^{(\mu\mu)} = -\kappa\cos(y-\mu), \qquad \ell^{(\mu\kappa)} =
\sin(y-\mu), \qquad \ell^{(\kappa\kappa)} = -A'(\kappa),\$\$ with
\\A(\kappa) = I_1(\kappa)/I_0(\kappa)\\. The pure-concentration entry is
**free of the data**: the log-density is \\\kappa\cos(y-\mu) -
\log\\2\pi I_0(\kappa)\\\\, linear in \\\kappa\\ apart from the
normalizing constant, so \\\kappa\\ appears twice only inside \\\log
I_0\\. It therefore equals its own expectation at every observation.

\\A'\\ comes from the Riccati recurrence \\A' = 1 - A/\kappa - A^2\\,
which
[`numericals7::bessel_i_ratio_derivs()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_derivs.html)
runs, so no second Bessel evaluation is needed. It is the variance of
\\\cos(Y-\mu)\\ and is positive, so the information is positive
definite.

## Arguments

- distrib:

  A `VonMises1Distrib` object, from
  [`vonmises1_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md).

- y:

  A numeric vector of angles in \\\[-\pi, \pi)\\.

- theta:

  A named list with components `mu` and `kappa`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `mu` must lie in \\(-\pi, \pi)\\ and `kappa` be strictly
  positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of three numeric vectors, `mu_mu`, `mu_kappa` and
`kappa_kappa`, each of length
`max(length(y), length(mu), length(kappa))`.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean
direction, \\\kappa \> 0\\ the concentration, \\I_m\\ the modified
Bessel function of the first kind of order \\m\\, and \\A(\kappa) =
I_1(\kappa)/I_0(\kappa)\\.

## See also

[`distrib_gradient.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.VonMises1Distrib.md)
for the score,
[`distrib_expected_hessian.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.VonMises1Distrib.md)
for the expectation of this quantity,
[`numericals7::bessel_i_ratio_derivs()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_derivs.html)
for \\A'\\, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- vonmises1_distrib()
y <- c(-1, 0, 0.5, 2)
th <- list(mu = 0.5, kappa = 2)
h <- distrib_hessian(d, y, th)
names(h)
#> [1] "mu_mu"       "mu_kappa"    "kappa_kappa"

# The pure-concentration entry is one number, repeated.
unique(h$kappa_kappa)
#> [1] -0.1642232

# And it is minus A', which the Riccati recurrence gives.
A <- numericals7::bessel_i_ratio(2)
c(riccati = 1 - A / 2 - A^2,
  supplied = numericals7::bessel_i_ratio_derivs(2)$d1)
#>   riccati  supplied 
#> 0.1642232 0.1642232 

# numDeriv on the summed log-density reproduces the summed matrix.
fn <- function(p)
  sum(distrib_pdf(d, y, list(mu = p[1], kappa = p[2]), log = TRUE))
H <- numDeriv::hessian(fn, c(0.5, 2))
rbind(numeric = c(H[1, 1], H[2, 2], H[1, 2]),
      closed = c(sum(h$mu_mu), sum(h$kappa_kappa), sum(h$mu_kappa)))
#>              [,1]       [,2]       [,3]
#> numeric -4.038114 -0.6568928 -0.4794255
#> closed  -4.038114 -0.6568928 -0.4794255
```
