# Render a Duration With a Unit Matched to Its Size

Formats a time in seconds using milliseconds below a second, seconds
below a minute, and minutes and seconds above, so that the figure a fit
prints stays readable whatever the sample size.

## Usage

``` r
fit_format_elapsed(sec)
```

## Arguments

- sec:

  A single non-negative number of seconds.

## Value

A character string.

## See also

[`print.distrib_fit`](https://statmodels7.github.io/distributions7/reference/print.distrib_fit.md)
