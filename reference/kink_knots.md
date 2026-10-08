# Where a Family's Kinks Sit on the Response Scale

Returns the points of \\\[\mathrm{lower}, \mathrm{upper}\]\\ at which
the log-density of a family is not smooth in its parameters, for one
parameter combination.
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
adds them to the knots of its quadrature, so that no panel straddles a
jump of the score.

## Usage

``` r
kink_knots(distrib, theta, lower, upper, ...)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- theta:

  A named list of parameter values, one value each.

- lower, upper:

  The finite range searched.

- ...:

  Unused.

## Value

A numeric vector, possibly empty.

## Details

The base method solves \\v(y, \theta) = 0\\ for the \\v\\ of the
family's
[`kink_decomposition()`](https://statmodels7.github.io/distributions7/reference/kink_decomposition.md).
It evaluates \\v\\ on a grid of the range, which is twice as dense on
its central half, and refines each sign change with
[`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html), so that a
\\v\\ with several roots (a Huber likelihood has two) gives each of
them. A family that declares no decomposition returns no point, even
when `params_smooth` records a non-smooth parameter, since the location
of its kink is then unknown.

A
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md)
family has a kink at \\\lvert y^\* \rvert\\ for every kink \\y^\*\\ of
its parent, its density being \\f(y) + f(-y)\\, and a
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
family has its parent's kinks.

## See also

[`kink_decomposition()`](https://statmodels7.github.io/distributions7/reference/kink_decomposition.md),
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md).

## Examples

``` r
distributions7:::kink_knots(laplace_distrib(), list(mu = 0.5, sigma = 2), -10, 10)
#> [1] 0.5
distributions7:::kink_knots(folded(laplace_distrib()),
                            list(mu = -0.5, sigma = 2), 0, 10)
#> [1] 0.5
```
