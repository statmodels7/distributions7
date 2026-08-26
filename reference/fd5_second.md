# A Five-Point Second Derivative

Returns \\\\-f(x-2h) + 16f(x-h) - 30f(x) + 16f(x+h) -
f(x+2h)\\/(12h^2)\\, the central stencil on five nodes for the second
derivative, with truncation error \\O(h^4)\\. Like
[`fd5_first()`](https://statmodels7.github.io/distributions7/reference/fd5_first.md)
it comes from
[`numericals7::fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.html)
at `accuracy = 4`.

Rounding is amplified by \\h^{-2}\\ here, one power more than in the
first derivative, so the attainable accuracy at the same step is one or
two digits lower.

## Usage

``` r
fd5_second(f, x, h)
```

## Arguments

- f:

  A function of one scalar, returning a numeric vector. It is called
  five times, at \\x\\, \\x \pm h\\ and \\x \pm 2h\\.

- x:

  A single number, the point to differentiate at.

- h:

  A single positive number, the step.

## Value

A numeric vector, of whatever length `f` returns.

## See also

[`fd5_first()`](https://statmodels7.github.io/distributions7/reference/fd5_first.md)
for the order below,
[`fd5_third()`](https://statmodels7.github.io/distributions7/reference/fd5_third.md)
for the order above, and
[`distrib_hessian.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.SkewTDistrib.md),
which uses this for the \\\nu\\ components.

## Examples

``` r
f <- function(x) exp(x) * sin(x)
c(stencil = distributions7:::fd5_second(f, 0.7, 1e-3),
  truth = 2 * exp(0.7) * cos(0.7))
#>  stencil    truth 
#> 3.080406 3.080406 
```
