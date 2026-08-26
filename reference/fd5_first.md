# A Five-Point First Derivative

Returns \\\\f(x-2h) - 8f(x-h) + 8f(x+h) - f(x+2h)\\/(12h)\\, the central
stencil on five nodes that is exact on polynomials up to degree four, so
its truncation error is \\O(h^4)\\. The weights come from
[`numericals7::fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.html)
at `accuracy = 4`, which builds them from the Vandermonde system on the
offsets \\-2, -1, 0, 1, 2\\.

## Usage

``` r
fd5_first(f, x, h)
```

## Arguments

- f:

  A function of one scalar, returning a numeric vector. It is called
  four times, at \\x \pm h\\ and \\x \pm 2h\\.

- x:

  A single number, the point to differentiate at.

- h:

  A single positive number, the step. For this family it comes from
  [`skewt_nu_step()`](https://statmodels7.github.io/distributions7/reference/skewt_nu_step.md).

## Value

A numeric vector, of whatever length `f` returns.

## Details

The three-point stencil is not accurate enough for the derivative in
\\\nu\\ of a fitted likelihood. Its truncation error is \\O(h^2)\\ per
observation and does **not** cancel when the observations are summed,
because it is a bias: on a sample of a few thousand it leaves the summed
score at about \\10^{-8}\\, an order of magnitude worse than the
five-point stencil, which reaches the level of rounding at the cost of
two more evaluations of `f`.

This is one stencil applied to an analytic quantity. Nothing in this
family differences a differenced value.

## See also

[`fd5_second()`](https://statmodels7.github.io/distributions7/reference/fd5_second.md),
[`fd5_third()`](https://statmodels7.github.io/distributions7/reference/fd5_third.md)
and
[`fd5_fourth()`](https://statmodels7.github.io/distributions7/reference/fd5_fourth.md)
for the other orders,
[`skewt_nu_step()`](https://statmodels7.github.io/distributions7/reference/skewt_nu_step.md)
for the step, and
[`numericals7::fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.html)
for the stencil library.

## Examples

``` r
f <- function(x) exp(x) * sin(x)
truth <- exp(0.7) * (sin(0.7) + cos(0.7))
c(stencil = distributions7:::fd5_first(f, 0.7, 1e-3), truth = truth)
#>  stencil    truth 
#> 2.837498 2.837498 

# The error falls as h^4 until rounding takes over.
vapply(c(0.1, 0.05, 0.025),
       function(h) abs(distributions7:::fd5_first(f, 0.7, h) - truth), 0)
#> [1] 3.784091e-05 2.364702e-06 1.477882e-07
```
