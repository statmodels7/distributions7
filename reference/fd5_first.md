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
leaves the summed score at about \\10^{-9}\\, which is above the default
gradient tolerance of
[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
and would make the fit run to its iteration limit and report failure at
a point that is in fact the maximum. The five-point stencil is
\\O(h^4)\\ and takes that to the level of rounding, at the cost of two
more evaluations.

This is one stencil applied to an analytic quantity, not a difference of
a difference.
