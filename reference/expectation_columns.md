# Aligned Parameter Columns for an Expectation

Shared preparation for the
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
methods: checks that the names in `...` do not collide with those of
`theta`, then expands every component to one aligned column per
parameter combination.

## Usage

``` r
expectation_columns(f_env_theta, dots)
```

## Arguments

- f_env_theta:

  The named list of parameters.

- dots:

  The list of further arguments destined for `f`.

## Value

A list with the theta columns `th`, the dot columns `dots` and the
number of combinations `n`.
