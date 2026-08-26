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

The identity is exact and needs no formula of its own, which is the
point: the location derivative is already written, often as a compiled
kernel, and the response derivative is the same number with a sign. It
also fixes the two orders below, where \\\partial\ell/\partial y =
-\partial\ell/ \partial\mu\\ and \\\partial^{2}\ell/\partial y^{2} =
\partial^{2}\ell/\partial\mu^{2}\\, and the tests check the new orders
against exactly that.

## Registered on

This body serves both third- and fourth-order response derivatives on
all fourteen location families, so `?distrib_deriv3_y.Gaussian1Distrib`
and its twenty-seven siblings open this page: `Gaussian1Distrib`,
`Gaussian2Distrib`, `Gaussian3Distrib`, `CauchyDistrib`,
`LogisticDistrib`, `LaplaceDistrib`, `Laplace2Distrib`, `EnetDistrib`,
`PseudoHuberDistrib`, `StudentT1Distrib`, `SkewNormal1Distrib`,
`SkewNormal2Distrib`, `SkewTDistrib`, `GumbelDistrib`.

## See also

[`distrib_deriv3_y.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.continuous_distrib.md)
for the stencil a family outside this list falls back to;
[`register_dy_k()`](https://statmodels7.github.io/distributions7/reference/register_dy_k.md),
the companion for families whose response is not a pure location;
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md),
whose components this reads.
