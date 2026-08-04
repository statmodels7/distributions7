# A Five-Point Third Derivative

Returns \\(-f(x-2h)/2 + f(x-h) - f(x+h) + f(x+2h)/2)/h^3\\, the
second-order central stencil for the third derivative. One stencil
applied to an analytic quantity, never a difference of differences.

## Usage

``` r
fd5_third(f, x, h)
```

## Arguments

- f:

  A function of one scalar, returning a numeric vector.

- x:

  The point to differentiate at.

- h:

  The step.

## Value

A numeric vector.
