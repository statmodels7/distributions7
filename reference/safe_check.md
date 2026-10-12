# Run a Check, Turning an Error Into a Failure

Evaluates a check expression and converts any error into a failed row
rather than letting it abort the report.

## Usage

``` r
safe_check(name, expr)
```

## Arguments

- name:

  The check's name, used for the row built on failure.

- expr:

  The expression to evaluate; normally returns a row from
  [`new_check()`](https://statmodels7.github.io/distributions7/reference/new_check.md).

## Value

The value of `expr`, or a failed row carrying the error message and the
attribute `from_error = TRUE`.

## Details

A distribution under validation may be defective, so an error raised by
a check is a result of that check. Without this conversion, the first
component that raised an error would end the run and the checks after it
would not be reported.

The failed row is marked with the attribute `from_error`, which
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
reads before assembling the table: a check whose computation raised an
error has failed, and it is not treated as a check whose statistic came
out `NaN` or `NA`, which is reported as not run.

## See also

[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md),
[`new_check()`](https://statmodels7.github.io/distributions7/reference/new_check.md)
