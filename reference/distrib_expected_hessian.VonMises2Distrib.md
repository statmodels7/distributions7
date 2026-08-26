# von Mises Expected Hessian in the Resultant Length

Returns the expectation of the observed Hessian under the model, in
closed form and with no quadrature or simulation. Since
\\\mathbb{E}\[\cos(Y-\mu)\] = A(\kappa)\\ and
\\\mathbb{E}\[\sin(Y-\mu)\] = 0\\, the term carrying the map's second
derivative drops out and \$\$\mathbb{E}\[\ell^{(\mu\mu)}\] = -\kappa
A(\kappa), \qquad \mathbb{E}\[\ell^{(\mu\rho)}\] = 0, \qquad
\mathbb{E}\[\ell^{(\rho\rho)}\] = -\dfrac{1}{A'(\kappa)}.\$\$

The last equality is the reparametrization identity: the information in
\\\kappa\\ is \\A'(\kappa)\\ and \\\kappa'(\rho) = 1/A'(\kappa)\\, so
\\A'(\kappa)\\\kappa'(\rho)\\^2 = 1/A'(\kappa)\\. The information in
\\\rho\\ is therefore the **reciprocal** of the information in
\\\kappa\\, which a one-to-one reparametrization of a single parameter
whose Jacobian is that reciprocal must give. The two parameters stay
orthogonal.

`approx` and `nsim` are ignored, the expectation being exact.

## Arguments

- distrib:

  A `VonMises2Distrib` object, from
  [`vonmises2_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises2_distrib.md).

- y:

  A numeric vector of angles. Only its length is read, the expectation
  not depending on the data; the values are ignored.

- theta:

  A named list with components `mu` and `rho`, each a numeric vector of
  length 1. `rho` must lie in \\(0, 1)\\.

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

A named list of three numeric vectors, `mu_mu`, `rho_rho` and `mu_rho`,
each of length `length(y)` and constant along it. `mu_rho` is exactly
zero.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean
direction, \\\rho \in (0,1)\\ the mean resultant length, \\\kappa\\ the
concentration and \\A(\kappa) = I_1(\kappa)/I_0(\kappa)\\.

## See also

[`distrib_hessian.VonMises2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.VonMises2Distrib.md)
for the quantity this is the expectation of,
[`distrib_expected_hessian.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.VonMises1Distrib.md)
for the same quantity in the concentration, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d2 <- vonmises2_distrib()
th <- list(mu = 0.5, rho = 0.7)
eh <- distrib_expected_hessian(d2, c(-1, 0, 0.5, 2), th)
vapply(eh, function(v) v[1], numeric(1))
#>     mu_mu   rho_rho    mu_rho 
#> -1.409539 -6.158821  0.000000 

# The reparametrization identity: the information in rho is the reciprocal
# of the information in kappa.
k <- numericals7::bessel_i_ratio_inverse(0.7)$kappa
Ap <- numericals7::bessel_i_ratio_derivs(k)$d1
c(supplied = eh$rho_rho[1], reciprocal = -1 / Ap)
#>   supplied reciprocal 
#>  -6.158821  -6.158821 

# Averaging the observed Hessian over draws reaches the same three numbers.
set.seed(1)
z <- distrib_rng(d2, 3e5, th)
vapply(distrib_hessian(d2, z, th), mean, numeric(1))
#>        mu_mu      rho_rho       mu_rho 
#> -1.410191065 -6.148592554 -0.004161727 
```
