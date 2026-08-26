# A Duration Rendered With a Unit Matched to Its Size

Formats a time in seconds for the line
[`print.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/print.distrib_fit.md)
shows, choosing the unit from the size: milliseconds below a second,
seconds below a minute, and minutes with seconds above. A fit of a few
hundred observations then reads `40 ms` where one of ten million reads
`2 min 07 s`, and neither prints a figure the reader has to count zeros
in.

## Usage

``` r
fit_format_elapsed(sec)
```

## Arguments

- sec:

  A single non-negative number of seconds. `NA`, a non-finite value and
  a zero-length vector all give `NA_character_`; the value is not
  otherwise validated, and a negative number is formatted as it stands.

## Value

A character string of length 1: `"0.4 ms"`, `"40 ms"`, `"1.5 s"`,
`"59.4 s"`, `"1 min 01 s"`. Seconds are rounded to one decimal below a
minute and to the nearest second above it.

## See also

[`print.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/print.distrib_fit.md),
the only caller;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
which accumulates the figure over every starting value and every
fallback.
