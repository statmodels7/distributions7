# The Response Derivative of a Location Family

Builds the order-`k` response derivative of a family whose response
enters only as \\y - \mu\\, as \\(-1)^{k}\\ times the pure derivative in
the location.

## Usage

``` r
loc_deriv_y_k(order)
```

## Arguments

- order:

  The derivative order, 3 or 4.

## Value

A function suitable for registering as an S7 method: it takes
`(distrib, y, theta, ...)` and returns a numeric vector of length
`length(y)`.

## Details

The identity is exact and needs no formula of its own. The location
derivative is already written, often as a compiled kernel, and the
response derivative equals it up to the sign \\(-1)^{k}\\. The same
identity holds at the lower orders, where \\\partial\ell/\partial y =
-\partial\ell/ \partial\mu\\ and \\\partial^{2}\ell/\partial y^{2} =
\partial^{2}\ell/\partial\mu^{2}\\.

## Registered on

This body serves both third- and fourth-order response derivatives on
the location families `Gaussian1Distrib`, `Gaussian2Distrib`,
`Gaussian3Distrib`, `CauchyDistrib`, `LogisticDistrib`,
`LaplaceDistrib`, `Laplace2Distrib`, `EnetDistrib`,
`PseudoHuberDistrib`, `StudentT1Distrib`, `SkewNormal1Distrib`,
`SkewTDistrib` and `GumbelDistrib`, and on `PseudoHuber2Distrib`, whose
page documents it. `SkewNormal2Distrib` has compiled kernels of its own
in the response; its two methods are listed here as aliases.

## See also

[`distrib_deriv3_y.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.continuous_distrib.md)
for the stencil a family outside this list falls back to;
[`register_dy_k()`](https://statmodels7.github.io/distributions7/reference/register_dy_k.md),
the companion for families whose response is not a pure location;
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md),
whose components this reads.
