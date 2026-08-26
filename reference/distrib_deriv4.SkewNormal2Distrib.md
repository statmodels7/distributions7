# Skew Normal Fourth Derivatives in the Centered Parametrization

Computes the fifteen fourth derivatives of the log-density in the
centered parameters, by the fourth-order partition sum of
[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md)
over the map of
[`sn_cp_to_dp()`](https://statmodels7.github.io/distributions7/reference/sn_cp_to_dp.md).
The map's fourth partial derivatives are the last entries of
[`md_skewnormal2()`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md)'s
table, so the order costs one more term of the same enumeration.

With `expected = TRUE` the parent's expected derivatives are carried
instead, and those are numerical.

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

  One of `"integrate"`, `"bartlett"`, `"mc"` or `"opg"`, read when
  `expected = TRUE`.

- nsim:

  A single positive integer, the Monte Carlo sample size used when
  `approx = "mc"`. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of fifteen numeric vectors, one per distinct fourth-order
component in the centered parameters.

## Errors

Signals an error when any element of `gamma1` is exactly zero.

## See also

[`distrib_deriv3.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.SkewNormal2Distrib.md)
for the order below,
[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md)
for the partition sum, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- skewnormal2_distrib()
y <- c(-1, 0.3, 1.7)
th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
length(distrib_deriv4(d, y, th))
#> [1] 15

# Against a central difference of the third order.
eps <- 1e-5
rbind(analytic = distrib_deriv4(d, y, th)$mu_mu_gamma1_gamma1,
      numeric = (distrib_deriv3(d, y, list(mu = 0, sigma = 1,
                                           gamma1 = 0.5 + eps))$mu_mu_gamma1 -
                 distrib_deriv3(d, y, list(mu = 0, sigma = 1,
                                           gamma1 = 0.5 - eps))$mu_mu_gamma1) /
                (2 * eps))
#>              [,1]     [,2]      [,3]
#> analytic -10.5455 2.616155 -1.060101
#> numeric  -10.5455 2.616155 -1.060101
```
