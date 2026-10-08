# The Roots of a Kink's Argument in a Range

Evaluates \\v(y, \theta)\\ of a
[`kink_spec()`](https://statmodels7.github.io/distributions7/reference/kink_spec.md)
on 129 equally spaced points of \\\[\mathrm{lower}, \mathrm{upper}\]\\
and 129 more on its central half, and refines each sign change with
[`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html). A grid point
at which \\v\\ is exactly zero is returned as it is.

## Usage

``` r
.kink_roots(kd, theta, lower, upper)
```

## Arguments

- kd:

  A
  [`kink_spec()`](https://statmodels7.github.io/distributions7/reference/kink_spec.md).

- theta:

  A named list of parameter values, one value each.

- lower, upper:

  The finite range searched.

## Value

A numeric vector of roots, possibly empty.

## Examples

``` r
d <- laplace_distrib()
distributions7:::.kink_roots(kink_decomposition(d), list(mu = 1, sigma = 1), -5, 5)
#> [1] 1
```
