# Where the Kink Sits on the Response Scale

Solves \\v(y, \theta) = 0\\ for `y`, which is where the log-density is
not smooth. The search runs over the support, cut to \\\[-1000, 1000\]\\
on a side whose bound is infinite, and returns `NA` when no sign change
is bracketed there.

## Usage

``` r
.kink_locate(distrib, kd, theta)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- kd:

  A
  [`kink_spec()`](https://statmodels7.github.io/distributions7/reference/kink_spec.md).

- theta:

  A named list of parameter values.

## Value

A single number, or `NA_real_`.

## Examples

``` r
d <- laplace_distrib()
distributions7:::.kink_locate(d, kink_decomposition(d), list(mu = 0.5, sigma = 2))
#> [1] 0.5
```
