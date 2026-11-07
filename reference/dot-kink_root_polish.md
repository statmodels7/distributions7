# Polish a Root of a Kink's Argument to Adjacent Doubles

Bisects the bracket \\\[a, b\]\\ on which \\v\\ changes sign, starting
from the root `r` that
[`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html) returned,
until the bracket is two adjacent doubles or a point with \\v = 0\\ is
met, and returns that point or the end of the bracket with the smaller
\\\|v\|\\. A root that is representable exactly (\\y = \mu\\ for \\v =
y - \mu\\) is then returned exactly, so that it coincides with a
quantile knot at the same point.

## Usage

``` r
.kink_root_polish(f, a, b, fa, fb, r)
```

## Arguments

- f:

  The function \\v\\ of one argument.

- a, b:

  The bracket.

- fa, fb:

  The values of `f` at `a` and `b`, of opposite signs.

- r:

  The root returned by
  [`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html).

## Value

A single number; `r` when `f` is not finite at a point visited.

## Examples

``` r
f <- function(x) x - 0.1
distributions7:::.kink_root_polish(f, 0, 1, -0.1, 0.9, 0.1 + 1e-15) == 0.1
#> [1] TRUE
```
