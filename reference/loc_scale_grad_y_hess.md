# Second-Order Mixed Derivatives of a Location-Scale Family

The
[`distrib_grad_y_hess`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
and
[`distrib_hess_y_hess`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
bodies shared by the families that are location-scale in both their
parameters.

## Usage

``` r
loc_scale_grad_y_hess(distrib, y, theta, scale = c("parameter", "link"), ...)

loc_scale_hess_y_hess(distrib, y, theta, scale = c("parameter", "link"), ...)
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

A named list keyed by parameter pair.

## See also

[`loc_scale_theta2_block`](https://statmodels7.github.io/distributions7/reference/loc_scale_theta2_block.md)
