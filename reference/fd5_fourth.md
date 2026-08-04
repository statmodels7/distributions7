# A Five-Point Fourth Derivative

Returns \\(f(x-2h) - 4f(x-h) + 6f(x) - 4f(x+h) + f(x+2h))/h^4\\, the
second-order central stencil for the fourth derivative. Rounding is
amplified by \\h^{-4}\\, so with the family's relative step this is
accurate to roughly four significant digits, which the pages that rely
on it state.

## Usage

``` r
fd5_fourth(f, x, h)
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
