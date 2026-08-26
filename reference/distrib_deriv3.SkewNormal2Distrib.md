# Skew Normal Third Derivatives in the Centered Parametrization

Computes the ten third derivatives of the log-density in the centered
parameters, by the third-order partition sum of
[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md)
over the map of
[`sn_cp_to_dp()`](https://statmodels7.github.io/distributions7/reference/sn_cp_to_dp.md).
The parent's third derivatives are closed form in a compiled kernel and
the map's partial derivatives are a written-out table in
[`md_skewnormal2()`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md),
so nothing here is a finite difference.

With `expected = TRUE` the parent's expected derivatives are carried
instead. Those are numerical, so `approx` and `nsim` are read.

## Arguments

- distrib:

  A `SkewNormal2Distrib` object, from
  [`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma` and `gamma1`. The skewness
  must not be exactly zero.

- expected:

  Logical of length 1. When `FALSE`, the default, the observed
  derivatives at `y` are returned.

- scale:

  Either `"parameter"`, the default, or `"link"`. The transformation is
  applied in the generic's body.

- approx:

  One of `"integrate"`, `"bartlett"`, `"mc"` or `"opg"`, the strategy
  the parent uses when `expected = TRUE`.

- nsim:

  A single positive integer, the Monte Carlo sample size used when
  `approx = "mc"`. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of ten numeric vectors, one per distinct third-order
component in the centered parameters, from `mu_mu_mu` to
`gamma1_gamma1_gamma1` as
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
names them.

## Errors

Signals an error when any element of `gamma1` is exactly zero.

## See also

[`distrib_hessian.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.SkewNormal2Distrib.md)
for the order below,
[`distrib_deriv4.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.SkewNormal2Distrib.md)
for the order above,
[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md)
for the partition sum, and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- skewnormal2_distrib()
y <- c(-1, 0.3, 1.7)
th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
d3 <- distrib_deriv3(d, y, th)
names(d3)
#>  [1] "mu_mu_mu"             "mu_mu_sigma"          "mu_mu_gamma1"        
#>  [4] "mu_sigma_sigma"       "mu_sigma_gamma1"      "mu_gamma1_gamma1"    
#>  [7] "sigma_sigma_sigma"    "sigma_sigma_gamma1"   "sigma_gamma1_gamma1" 
#> [10] "gamma1_gamma1_gamma1"

# Against a central difference of the analytic Hessian.
eps <- 1e-5
rbind(analytic = d3$mu_mu_gamma1,
      numeric = (distrib_hessian(d, y, list(mu = 0, sigma = 1,
                                            gamma1 = 0.5 + eps))$mu_mu -
                 distrib_hessian(d, y, list(mu = 0, sigma = 1,
                                            gamma1 = 0.5 - eps))$mu_mu) /
                (2 * eps))
#>               [,1]     [,2]      [,3]
#> analytic -3.132477 1.233707 0.3509893
#> numeric  -3.132477 1.233707 0.3509893
```
