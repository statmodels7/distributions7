# Colors and Line Types for Overlaid Curves

The visual keys distinguishing several settings on one panel: a color
and a line type per curve, both cycled, so the curves are told apart in
color and in a printed copy that has none.

## Usage

``` r
plot_keys(k, dots = list())
```

## Arguments

- k:

  The number of curves.

- dots:

  The caller's `...`; a `col` or `lty` given there wins and is recycled
  over the curves.

## Value

A list with `col` and `lty`, each of length `k`.
