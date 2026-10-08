# A Family at Selected Observations

Returns `distrib` with every constant that varies by observation (a
binomial's or a beta-binomial's `size`), its parents' included, taken at
the observations `idx`, so that the family can be evaluated at points
that belong to those observations.

## Usage

``` r
distrib_at_rows(distrib, idx, n)
```

## Arguments

- distrib:

  A univariate family.

- idx:

  Integer indices of observations, one per point.

- n:

  The number of observations the constants are recycled to.

## Value

The family, with its varying constants of length `length(idx)`.
