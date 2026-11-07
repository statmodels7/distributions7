# The Constants of a Family That Vary by Observation

Collects the constants of `distrib` and of its parent families that vary
by observation (a binomial's `size`), each recycled to `n`.

## Usage

``` r
distrib_const_columns(distrib, n)
```

## Arguments

- distrib:

  A univariate family.

- n:

  The number of observations.

## Value

A list of numeric vectors of length `n`, empty when the family has no
such constant.

## Examples

``` r
distributions7:::distrib_const_columns(binomial_distrib(size = c(5, 10)), 4)
#> [[1]]
#> [1]  5 10  5 10
#> 
```
