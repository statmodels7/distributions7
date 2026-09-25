# Pseudo-Huber Expected Hessian and Its Derivatives

Returns the expectation of the observed Hessian under the model, and
through
[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)
and
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md)
its first and second derivatives in the parameters.

## Arguments

- distrib:

  A `PseudoHuberDistrib` object, from
  [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).

- y:

  A numeric vector of observations. Its length sets the length of each
  returned component; the values themselves are not read.

- theta:

  A named list with components `mu`, `sigma` and `nu`, each a numeric
  vector of length 1 or of the length of `y`. `sigma` and `nu` must be
  strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html).

- approx, nsim:

  Ignored.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  The thread count passed to the family's kernels.

## Value

A named list of numeric vectors of length `length(y)`: for
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
the six components `mu_mu`, `sigma_sigma`, `nu_nu`, `mu_sigma`, `mu_nu`
and `sigma_nu`; for the two derivatives the components keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
and
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md).

## Details

The family is a location-scale family, so every component equals its
value at \\\mu = 0\\, \\\sigma = 1\\ times \\\sigma^{-k}\\, with \\k\\
the number of indices on \\\mu\\ or \\\sigma\\. The value at the
standard location and scale is a function of \\\nu\\ alone and is an
integral over \\z\\ of the analytic observed derivatives against the
density, taken once for each distinct \\\nu\\ by the exp-sinh rule of
[`loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_expected.md).
No component has an elementary closed form: the integrands carry
\\(\nu + z^2)^{-1/2}\\, which leads to Bickley functions rather than to
the Bessel functions of the normalizing constant.

The law is symmetric about \\\mu\\, so every component carrying \\\mu\\
an odd number of times vanishes. The rule's nodes are symmetric about
zero, so the two halves of such an integral cancel and the entry comes
back as zero.

`approx` and `nsim` are accepted for the generic's sake and ignored.

## See also

[`loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_expected.md)
for the construction,
[`distrib_hessian.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.PseudoHuberDistrib.md)
for the quantity this is the expectation of, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- pseudohuber_distrib()
y <- c(-2.5, 0.3, 1.8)
th <- list(mu = 0.4, sigma = 1.2, nu = 2)
eh <- distrib_expected_hessian(d, y, th)
vapply(eh, function(v) v[1], numeric(1))
#>         mu_mu   sigma_sigma         nu_nu      mu_sigma         mu_nu 
#> -2.620828e-01 -9.177697e-01 -5.412069e-03 -1.616963e-20  1.769121e-21 
#>      sigma_nu 
#> -6.699758e-02 

# The entries odd in the residual vanish by symmetry.
c(eh$mu_sigma[1], eh$mu_nu[1])
#> [1] -1.616963e-20  1.769121e-21

# The information does not depend on the observations, and scales with
# sigma^-2 in the location and scale block.
e2 <- distrib_expected_hessian(d, y, list(mu = 0.4, sigma = 2.4, nu = 2))
eh$mu_mu[1] / e2$mu_mu[1]
#> [1] 4
```
