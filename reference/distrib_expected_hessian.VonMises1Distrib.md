# von Mises Expected Hessian

Returns the expectation of the observed Hessian under the model, in
closed form and with no quadrature or simulation. Under the model
\\\mathbb{E}\[\cos(Y-\mu)\] = A(\kappa)\\ and
\\\mathbb{E}\[\sin(Y-\mu)\] = 0\\ by symmetry, so
\$\$\mathbb{E}\[\ell^{(\mu\mu)}\] = -\kappa A(\kappa), \qquad
\mathbb{E}\[\ell^{(\mu\kappa)}\] = 0, \qquad
\mathbb{E}\[\ell^{(\kappa\kappa)}\] = -A'(\kappa),\$\$ with \\A(\kappa)
= I_1(\kappa)/I_0(\kappa)\\.

The zero off-diagonal makes the direction and the concentration
**orthogonal**: the expected information is diagonal, their estimates
are asymptotically independent, and Fisher scoring updates the two
independently. `approx` and `nsim` are ignored, the expectation being
exact.

## Arguments

- distrib:

  A `VonMises1Distrib` object, from
  [`vonmises1_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md).

- y:

  A numeric vector of angles. Only its length is read, the expectation
  not depending on the data; the values are ignored.

- theta:

  A named list with components `mu` and `kappa`, each a numeric vector
  of length 1. `kappa` must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored here, the expectation being exact. Accepted so that the
  signature matches the generic's, where it selects between
  `"bartlett"`, `"integrate"`, `"mc"` and `"opg"`.

- nsim:

  Ignored here, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of three numeric vectors, `mu_mu`, `mu_kappa` and
`kappa_kappa`, each of length `length(y)` and constant along it.
`mu_kappa` is exactly zero.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean
direction, \\\kappa \> 0\\ the concentration, and \\A(\kappa) =
I_1(\kappa)/I_0(\kappa)\\ the mean resultant length.

## See also

[`distrib_hessian.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.VonMises1Distrib.md)
for the quantity this is the expectation of,
[`distrib_gradient.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.VonMises1Distrib.md)
for the score, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- vonmises1_distrib()
th <- list(mu = 0.5, kappa = 2)
eh <- distrib_expected_hessian(d, c(-1, 0, 0.5, 2), th)
vapply(eh, function(v) v[1], numeric(1))
#>       mu_mu    mu_kappa kappa_kappa 
#>  -1.3955493   0.0000000  -0.1642232 

# The closed forms, written out; the off-diagonal is exactly zero, so the
# two parameters are orthogonal.
A <- numericals7::bessel_i_ratio(2)
c(mu_mu = -2 * A, mu_kappa = 0,
  kappa_kappa = -numericals7::bessel_i_ratio_derivs(2)$d1)
#>       mu_mu    mu_kappa kappa_kappa 
#>  -1.3955493   0.0000000  -0.1642232 

# Averaging the observed Hessian over draws reaches the same three numbers.
set.seed(1)
z <- distrib_rng(d, 3e5, th)
vapply(distrib_hessian(d, z, th), mean, numeric(1))
#>         mu_mu      mu_kappa   kappa_kappa 
#> -1.3964000867 -0.0006979838 -0.1642231977 

# The strategy argument is inert, the expectation being exact.
identical(eh, distrib_expected_hessian(d, c(-1, 0, 0.5, 2), th,
                                       approx = "mc", nsim = 50))
#> [1] TRUE
```
