# A Five-Point Third Derivative

Returns \\\\-f(x-2h)/2 + f(x-h) - f(x+h) + f(x+2h)/2\\/h^3\\, the
central stencil on five nodes for the third derivative. Five nodes is
the smallest number that carries a third derivative at all, so the
accuracy here is \\O(h^2)\\ where the first two orders get \\O(h^4)\\
from the same offsets.

One stencil applied to an analytic quantity, never a difference of
differences.

## Usage

``` r
fd5_third(f, x, h)
```

## Arguments

- f:

  A function of one scalar, returning a numeric vector. It is called
  four times, at \\x \pm h\\ and \\x \pm 2h\\.

- x:

  A single number, the point to differentiate at.

- h:

  A single positive number, the step.

## Value

A numeric vector, of whatever length `f` returns.

## See also

[`fd5_second()`](https://statmodels7.github.io/distributions7/reference/fd5_second.md)
for the order below,
[`fd5_fourth()`](https://statmodels7.github.io/distributions7/reference/fd5_fourth.md)
for the order above, and
[`distrib_deriv3.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.SkewTDistrib.md),
which uses this.

## Examples

``` r
f <- function(x) exp(x) * sin(x)
c(stencil = distributions7:::fd5_third(f, 0.7, 1e-2),
  truth = 2 * exp(0.7) * (cos(0.7) - sin(0.7)))
#>   stencil     truth 
#> 0.4855321 0.4858158 
```
