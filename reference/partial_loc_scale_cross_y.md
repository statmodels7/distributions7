# Mixed Derivatives When Only Two Parameters Are Location-Scale

The location and scale components in closed form and the remaining shape
components by one central difference of
[`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md).

## Usage

``` r
partial_loc_scale_cross_y(
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

[`loc_scale_cross_block`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross_block.md)
