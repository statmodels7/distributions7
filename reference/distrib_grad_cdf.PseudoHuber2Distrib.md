# Pseudo-Huber Distribution-Function Derivatives, Standard-Deviation Parametrization

The derivatives of the distribution function in the parameters, by the
location-scale identities of
[`partial_loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_grad_cdf.md),
[`partial_loc_scale_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_hess_cdf.md)
and, at orders three and four,
[`partial_loc_scale_deriv_cdf_k()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_deriv_cdf_k.md):
the location and scale components are closed in the density, and the
components in \\\nu\\ are differences of the quadrature.

## Arguments

- distrib:

  A `PseudoHuber2Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu`, `sigma` and `nu`.

- lower.tail, log:

  As in
  [`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md).

- ...:

  Unused.

## Value

A named list, one component per parameter (first order) or per pair
(second order).

## See also

[`distrib_cdf.PseudoHuber2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.PseudoHuber2Distrib.md).

## Examples

``` r
d <- pseudohuber2_distrib()
th <- list(mu = 0.4, sigma = 1.5, nu = 2)
all.equal(distrib_grad_cdf(d, 1, th, log = FALSE)$mu,
          -distrib_pdf(d, 1, th))
#> [1] TRUE
```
