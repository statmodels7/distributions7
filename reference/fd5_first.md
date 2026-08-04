# A Five-Point First Derivative

Returns \\(f(x-2h) - 8f(x-h) + 8f(x+h) - f(x+2h))/(12h)\\, the
fourth-order central stencil.

## Usage

``` r
fd5_first(f, x, h)
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

## Details

The three-point stencil is not accurate enough for the derivative in
\\\nu\\ of a fitted likelihood. Its truncation error is \\O(h^2)\\ per
observation and does **not** cancel when the observations are summed,
because it is a bias rather than noise: on a sample of a few thousand it
leaves the summed score at about \\10^{-8}\\, an order of magnitude
worse than the five-point stencil, which is \\O(h^4)\\ and reaches the
level of rounding at the cost of two more evaluations.

This is one stencil applied to an analytic quantity, not a difference of
a difference.
