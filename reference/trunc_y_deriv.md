# Response Derivative of a Truncated Distribution

Evaluates a response derivative of the parent inside the interval.

## Usage

``` r
trunc_y_deriv(distrib, y, theta, fun)
```

## Arguments

- distrib:

  A truncated distribution object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- fun:

  The parent's response-derivative function.

## Value

A numeric vector.

## Details

The normalising constant does not depend on \\y\\, so inside the support
the parent's derivative is exact and nothing needs correcting.
