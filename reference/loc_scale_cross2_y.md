# Second-Response Mixed Derivatives of a Location-Scale Family

The
[`distrib_cross2_y`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
body shared by the families that are location-scale in both their
parameters, and the partial form for those with a shape parameter beyond
the two.

## Usage

``` r
loc_scale_cross2_y(distrib, y, theta, scale = c("parameter", "link"), ...)

partial_loc_scale_cross2_y(
  distrib,
  y,
  theta,
  scale = c("parameter", "link"),
  ...
)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- scale:

  Handled by the generic before dispatch.

- ...:

  Unused.

## Value

A named list with one numeric vector per parameter.

## See also

[`loc_scale_cross2_block`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross2_block.md)
