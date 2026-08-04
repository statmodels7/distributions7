# A Five-Point Second Derivative

Returns \\(-f(x-2h) + 16f(x-h) - 30f(x) + 16f(x+h) - f(x+2h))/(12h^2)\\,
the fourth-order central stencil for the second derivative.

## Usage

``` r
fd5_second(f, x, h)
```

## Arguments

- f:

  A function of one scalar, returning a numeric vector.

- x:

  The point to differentiate at.

- h:

  The step.

## Value

A numeric vector, of whatever length `f` returns.
