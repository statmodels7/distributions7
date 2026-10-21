# Number of Observations Over Which a Family's Constants Vary

Returns the length of the longest constant of `distrib` that varies by
observation (a binomial's `size`), including those of its parent
families, and 1 when the family has none.

## Usage

``` r
distrib_const_rows(distrib)
```

## Arguments

- distrib:

  A univariate family.

## Value

A single positive integer.
