# Set Aside the Rows of a Validation That Have No Verdict

Binds the rows of
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
into a table, leaving out each row that compared nothing and recording
it with its reason in the attribute `"skipped"`: a row marked as a
numerical fallback, and a row whose statistic is `NaN` or `NA` although
its computation did not raise an error. A row that raised keeps its
failure, and so does a row named in `on_values`, whose criterion is
defined on the values themselves.

## Usage

``` r
drop_no_verdict(res, on_values)
```

## Arguments

- res:

  A list of one-row data frames from
  [`new_check()`](https://statmodels7.github.io/distributions7/reference/new_check.md)
  and
  [`safe_check()`](https://statmodels7.github.io/distributions7/reference/safe_check.md).

- on_values:

  A character vector of check names that are kept whatever their
  statistic.

## Value

A data frame with columns `check`, `status`, `statistic` and `detail`,
carrying the attribute `"skipped"` (a data frame with columns `check`
and `reason`) when some row was left out.

## See also

[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md),
the caller.
