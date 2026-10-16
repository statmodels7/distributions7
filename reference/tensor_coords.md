# The Coordinates of a Tensor Stencil over the Parameters

Returns the coordinates on which
[`tensor_derivatives()`](https://statmodels7.github.io/distributions7/reference/tensor_derivatives.md),
[`mixed_tensor_derivatives()`](https://statmodels7.github.io/distributions7/reference/mixed_tensor_derivatives.md)
and
[`numerical_deriv5()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv5.md)
difference: the parameters themselves on the parameter scale, or the
linear predictors \\\eta = g(\theta)\\ on the link scale, with the
bounds of each and the map back to a `theta` that the family's methods
accept.

## Usage

``` r
tensor_coords(distrib, theta, scale = c("parameter", "link"))
```

## Arguments

- distrib:

  An object inheriting from `distrib`.

- theta:

  A named list of parameters, aligned by the generic.

- scale:

  `"parameter"` or `"link"`.

## Value

A list with `x`, the named list of coordinates; `bounds`, a list of
`c(lower, upper)` per coordinate; and `to_theta`, a function from a list
of coordinates to `theta`.
