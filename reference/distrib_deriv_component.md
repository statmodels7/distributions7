# One Component of the Parent's Derivative

Fetches a single component of the parent's derivative of a given order,
by canonical key.

## Usage

``` r
distrib_deriv_component(parent, y, theta, idx, params, order)
```

## Arguments

- parent:

  The parent distribution.

- y:

  A numeric vector of observations.

- theta:

  A named list of the parent's parameters.

- idx:

  A character vector of parameter names, the multi-index.

- params:

  The parent's parameter names, in declaration order.

- order:

  The derivative order.

## Value

A numeric vector.
