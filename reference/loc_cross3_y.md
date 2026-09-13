# The Mixed Third-Response Derivative of a Location Family

Reads \\\partial^4\ell/\partial y^3\\\partial\theta_i\\ as
\\-\partial^4\ell/\partial\mu^3\\\partial\theta_i\\ from
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md),
for a family whose response enters only as \\y - \mu\\.

## Usage

``` r
loc_cross3_y(distrib, y, theta, scale = c("parameter", "link"), ...)
```

## Arguments

- distrib:

  A location-family distribution object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- scale:

  One of `"parameter"` or `"link"`, applied by the generic after
  dispatch.

- ...:

  Unused.

## Value

A named list with one numeric vector per parameter, each as long as `y`.

## Details

The identity is exact: each derivative in the response is minus the
derivative in the location, and a derivative in any parameter commutes
with both. The component for \\\theta_i = \mu\\ is therefore
\\-\partial^4\ell/\partial\mu^4\\. The keys are generated from the
parameter indices, as
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
generates them, rather than written.

## Registered on

The fourteen families
[`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
serves through
[`loc_deriv_y_k()`](https://statmodels7.github.io/distributions7/reference/loc_deriv_y_k.md):
`Gaussian1Distrib`, `Gaussian2Distrib`, `Gaussian3Distrib`,
`CauchyDistrib`, `LogisticDistrib`, `LaplaceDistrib`, `Laplace2Distrib`,
`EnetDistrib`, `PseudoHuberDistrib`, `StudentT1Distrib`,
`SkewNormal1Distrib`, `SkewNormal2Distrib`, `SkewTDistrib`,
`GumbelDistrib`.

## See also

[`distrib_cross3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross3_y.md),
[`loc_deriv_y_k()`](https://statmodels7.github.io/distributions7/reference/loc_deriv_y_k.md).
