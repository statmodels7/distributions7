# Second-Order Mixed Derivatives When Only Two Parameters Are Location-Scale

The location and scale pairs in closed form and every pair touching a
shape parameter by one central difference of the analytic first-order
component.

## Usage

``` r
partial_loc_scale_grad_y_hess(
  distrib,
  y,
  theta,
  scale = c("parameter", "link"),
  ...
)

partial_loc_scale_hess_y_hess(
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

A named list keyed by parameter pair.

## See also

[`loc_scale_theta2_block`](https://statmodels7.github.io/distributions7/reference/loc_scale_theta2_block.md)
