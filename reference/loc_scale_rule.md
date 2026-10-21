# The Quadrature Rule of the Location-Scale Expectations

Returns the nodes and weights of an exp-sinh rule on each half-line, for
integrals over the whole real line of a density at location 0 and scale
1.

## Usage

``` r
loc_scale_rule()
```

## Value

A list with the nodes `x` and the weights `w`.

## Details

On \\(0,\infty)\\ the substitution \\z = \exp(\tfrac{\pi}{2}\sinh t)\\
turns the integral into one over the whole \\t\\ line whose integrand
decays double-exponentially at both ends, which the trapezoidal rule
then integrates with an error falling like \\\exp(-c/h)\\; the negative
half is its mirror image. The step is \\h = 1/32\\ and the nodes run
from \\10^{-16}\\ to \\10^{60}\\ in \\\|z\|\\, about 580 in all. Past
the upper end a tail decaying as \\\|z\|^{-1-\nu}\\ leaves a relative
\\10^{-60\nu}/\nu\\, negligible for \\\nu \ge 0.3\\.

## See also

[`loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_expected.md)
