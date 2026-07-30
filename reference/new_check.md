# Record One Check Result

Builds the single-row data frame that
[`check_distrib`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
accumulates into its report.

## Usage

``` r
new_check(name, ok, stat, detail = NA_character_)
```

## Arguments

- name:

  The check's name, as it appears in the report.

- ok:

  Logical; whether the check passed.

- stat:

  The numeric statistic the check produced, typically a maximum
  discrepancy.

- detail:

  An optional message, used to carry the reason for a failure.

## Value

A one-row data frame with columns `check`, `status`, `statistic` and
`detail`.

## See also

[`check_distrib`](https://statmodels7.github.io/distributions7/reference/check_distrib.md),
[`safe_check`](https://statmodels7.github.io/distributions7/reference/safe_check.md)
