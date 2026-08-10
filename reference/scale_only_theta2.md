# The Scale Component of a Family With No Location

\\\partial^2\ell^{(y)}/\partial\sigma^2\\ and
\\\partial^2\ell^{(yy)}/\partial\sigma^2\\ for a family whose response
enters only through \\z = y/\sigma\\.

## Usage

``` r
scale_only_theta2(distrib, y, theta, order = 1L, at = 1L)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- order:

  `1` for the block of
  [`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
  `2` for that of
  [`distrib_hess_y`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md).

- at:

  The index of the scale parameter.

## Value

A numeric vector.

## Details

The derivation of the scale components of
[`loc_scale_theta2_block`](https://statmodels7.github.io/distributions7/reference/loc_scale_theta2_block.md)
never used the location, only that \\\sigma\\ is a scale, so the
formulas hold unchanged with \\z = y/\sigma\\. Any other parameter is a
shape and is not covered.

## See also

[`loc_scale_theta2_block`](https://statmodels7.github.io/distributions7/reference/loc_scale_theta2_block.md)
