# The Beta-Binomial's Cumulative Mass Over Its Support

Returns the cumulative mass over \\\\0, \dots, n\\\\: a vector where
both parameters are scalar, and otherwise a matrix with one row per
observation, in which each observation's parameters are crossed with the
support. With `upper = TRUE` the entry at \\k\\ is \\P(Y \> k)\\, the
mass above \\k\\ accumulated from the top of the support, in place of
\\P(Y \le k)\\.

## Usage

``` r
betabinom1_cum(n, theta, len, upper = FALSE)
```

## Arguments

- n:

  The size.

- theta:

  A named list with the two parameters.

- len:

  The length of the quantity at which the caller evaluates the
  cumulative mass.

- upper:

  Logical of length 1. When `TRUE` the upper tail is returned. Defaults
  to `FALSE`.

## Value

A numeric vector of length `n + 1`, or a matrix with `n + 1` columns.

## Details

The kernel reads its parameters element by element, so the support and a
vector of parameters are not passed to it together: each observation's
parameters are crossed with the support, and the sum is taken along each
row.
