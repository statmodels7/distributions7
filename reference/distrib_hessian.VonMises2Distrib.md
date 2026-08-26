# von Mises Observed Hessian in the Resultant Length

Computes the three distinct second derivatives of the log-density in
\\\mu\\ and \\\rho\\, one value per observation, in closed form. The
concentration parametrization's second derivatives are carried through
the one-variable chain rule, \$\$\ell^{(\rho\rho)} =
\ell^{(\kappa\kappa)}\\\kappa'(\rho)\\^2 +
\ell^{(\kappa)}\kappa''(\rho), \qquad \ell^{(\mu\rho)} =
\ell^{(\mu\kappa)}\kappa'(\rho),\$\$ with \\\ell^{(\kappa\kappa)} =
-A'(\kappa)\\, \\\ell^{(\mu\kappa)} = \sin(y-\mu)\\ and
\\\ell^{(\mu\mu)} = -\kappa\cos(y-\mu)\\ unchanged.

Unlike in the concentration parametrization, the pure second derivative
is **not** free of the data: the term in \\\kappa''\\ carries
\\\cos(y-\mu)\\, which the map's curvature brings in.

## Arguments

- distrib:

  A `VonMises2Distrib` object, from
  [`vonmises2_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises2_distrib.md).

- y:

  A numeric vector of angles in \\\[-\pi, \pi)\\.

- theta:

  A named list with components `mu` and `rho`, each a numeric vector of
  length 1 or of the length of `y`. `mu` must lie in \\(-\pi, \pi)\\ and
  `rho` in \\(0, 1)\\.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of three numeric vectors, `mu_mu`, `rho_rho` and `mu_rho`,
each of length `max(length(y), length(mu), length(rho))`.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean
direction, \\\rho \in (0,1)\\ the mean resultant length, \\\kappa\\ the
concentration and \\A(\kappa) = I_1(\kappa)/I_0(\kappa)\\.

## See also

[`distrib_gradient.VonMises2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.VonMises2Distrib.md)
for the score,
[`distrib_expected_hessian.VonMises2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.VonMises2Distrib.md)
for the expectation of this quantity,
[`distrib_hessian.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.VonMises1Distrib.md)
for the same quantity in the concentration, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d2 <- vonmises2_distrib()
y <- c(-1, 0, 0.5, 2)
th <- list(mu = 0.5, rho = 0.7)
h <- distrib_hessian(d2, y, th)
names(h)
#> [1] "mu_mu"   "rho_rho" "mu_rho" 

# numDeriv on the summed log-density reproduces the summed matrix.
fn <- function(v)
  sum(distrib_pdf(d2, y, list(mu = v[1], rho = v[2]), log = TRUE))
H <- numDeriv::hessian(fn, c(0.5, 0.7))
rbind(numeric = c(H[1, 1], H[2, 2], H[1, 2]),
      closed = c(sum(h$mu_mu), sum(h$rho_rho), sum(h$mu_rho)))
#>              [,1]      [,2]      [,3]
#> numeric -4.065629 -49.32113 -2.952696
#> closed  -4.065629 -49.32113 -2.952696

# The pure second derivative varies with the data here, where in the
# concentration parametrization it does not.
h$rho_rho
#> [1] -26.0500105  -0.5453825   3.3242700 -26.0500105
```
