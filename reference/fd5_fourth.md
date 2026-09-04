# A Five-Point Fourth Derivative

Returns \\\\f(x-2h) - 4f(x-h) + 6f(x) - 4f(x+h) + f(x+2h)\\/h^4\\, the
central stencil on five nodes for the fourth derivative, accurate to
\\O(h^2)\\.

Rounding is amplified by \\h^{-4}\\, which sets what this can deliver:
measured on \\e^x\sin x\\ at \\x = 0.7\\ with \\h = 10^{-2}\\ it returns
\\-5.18939\\ against a true \\-5.18918\\, four significant digits. On a
component whose value is itself small the surviving digits are fewer,
and
[`distrib_deriv4.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.SkewTDistrib.md)
says so where a reader meets one.

## Usage

``` r
fd5_fourth(f, x, h)
```

## Arguments

- f:

  A function of one scalar, returning a numeric vector. It is called
  five times, at \\x\\, \\x \pm h\\ and \\x \pm 2h\\.

- x:

  A single number, the point to differentiate at.

- h:

  A single positive number, the step. A step too small is worse than one
  too large here, the rounding growing four times faster than the
  truncation falls.

## Value

A numeric vector, of whatever length `f` returns.

## See also

[`fd5_third()`](https://statmodels7.github.io/distributions7/reference/fd5_third.md)
for the order below and
[`distrib_deriv4.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.SkewTDistrib.md)
for the method that uses this.

## Examples

``` r
f <- function(x) exp(x) * sin(x)
truth <- -4 * exp(0.7) * sin(0.7)
c(stencil = distributions7:::fd5_fourth(f, 0.7, 1e-2), truth = truth)
#>   stencil     truth 
#> -5.189386 -5.189180 

# Too small a step is worse than too large: rounding grows as h^-4.
vapply(c(1e-1, 1e-2, 1e-3, 1e-4),
       function(h) abs(distributions7:::fd5_fourth(f, 0.7, h) - truth), 0)
#> [1] 0.0205100667 0.0002052736 0.0004421196 3.6926037495
```
