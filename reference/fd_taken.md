# Central Differences on the Steps Actually Taken

`fd_first_taken()` and `fd_second_taken()` form the central differences
of first and second order, dividing by the distances their evaluation
points actually lie at rather than by the nominal step.

## Usage

``` r
fd_first_taken(f, x, h)

fd_second_taken(f, x, h, f0 = f(x))
```

## Arguments

- f:

  A vectorized function of the points.

- x:

  The points, a numeric vector.

- h:

  The steps, as long as `x` or of length one.

- f0:

  `f(x)`, when the caller already has it.

## Value

A numeric vector as long as `x`.

## Details

The point \\x + h\\ is a double, so the step it lies at is \\h\_+ = (x +
h) - x\\, a subtraction that is exact, and the step on the other side is
\\h\_- = x - (x - h)\\. Both may differ from \\h\\, and from each other,
by the rounding of the evaluation points, at most one unit of the last
place. With them \$\$f'(x) \approx \frac{f(x+h) - f(x-h)}{h\_+ + h\_-},
\qquad f''(x) \approx \frac{2\\h\_- f(x+h) - (h\_+ + h\_-) f(x) + h\_+
f(x-h)\\} {h\_+ h\_- (h\_+ + h\_-)},\$\$ which are the uniform formulas
wherever \\h\_+ = h\_- = h\\. The correction is of the order of one unit
of the last place over the step: negligible for a step of many units,
and the whole of the accuracy for a step of a few, which is what a step
within about \\10^{-10}\\ of a non-zero bound is.

## See also

[`fd_stable_quotient()`](https://statmodels7.github.io/distributions7/reference/fd_stable_quotient.md),
whose callers pass these.

## Examples

``` r
# 1e-10 below the bound of log(1 - x), a step of 6e-16 is 5.4 units of the
# last place, and the nominal quotient carries the rounding of its points
x <- 1 - 1e-10
h <- 6e-16
f <- function(v) log1p(-v)
exact <- -1 / (1 - x)
c(nominal = (f(x + h) - f(x - h)) / (2 * h) / exact - 1,
  taken = distributions7:::fd_first_taken(f, x, h) / exact - 1)
#>       nominal         taken 
#> -7.481415e-02 -1.396359e-10 
```
