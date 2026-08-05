# Recycle a Constant Moment to the Length of the Parameters

Returns `value` repeated to the length the parameters imply, which is
what a moment that does not depend on them still has to be.

## Usage

``` r
moment_const(theta, k, value)
```

## Arguments

- theta:

  An aligned named list of parameters.

- k:

  The number of parameters to read the length from.

- value:

  The constant.

## Value

A numeric vector.
