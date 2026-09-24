# The Beta-Binomial's Cumulative Mass Over Its Support

The cumulative mass over \\\\0, \dots, n\\\\: a vector where both
parameters are scalar, and otherwise a matrix with one row per
observation, each observation's parameters crossed with the support.

## Usage

``` r
betabinom1_cum(n, theta, len)
```

## Arguments

- n:

  The size.

- theta:

  A named list with the two parameters.

- len:

  The length of the quantity the caller evaluates at.

## Value

A numeric vector of length `n + 1`, or a matrix with `n + 1` columns.

## Details

The kernel reads its parameters per element, so handing it the support
against a vector of parameters recycled nothing and read past the end of
the shorter vector; that is what the matrix branch replaced.
