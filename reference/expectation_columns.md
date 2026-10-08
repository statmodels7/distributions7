# Aligned Parameter Columns for an Expectation

Shared preparation for the
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
methods: checks that the names in `...` do not collide with those of
`theta`, then expands every component to one aligned column per
parameter combination. A constant of the family that varies by
observation (a binomial's `size`) counts as a column too, so that a
scalar `theta` gives one combination per observation; the methods take
the constant at each combination through
[`distrib_at_rows()`](https://statmodels7.github.io/distributions7/reference/distrib_at_rows.md).

## Usage

``` r
expectation_columns(distrib, f_env_theta, dots)
```

## Arguments

- distrib:

  The family.

- f_env_theta:

  The named list of parameters.

- dots:

  The list of further arguments destined for `f`.

## Value

A list with the theta columns `th`, the dot columns `dots`, the number
of combinations `n` and `rows`, which is `TRUE` when the family carries
a constant that varies by observation.
