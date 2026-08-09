# Colors, Line Types and Symbols for Overlaid Settings

The visual keys distinguishing several settings on one panel: a color, a
line type and a plotting symbol per setting, all cycled, so the settings
are told apart in color and in a printed copy that has none.

## Usage

``` r
plot_keys(k, dots = list())
```

## Arguments

- k:

  The number of settings.

- dots:

  The caller's `...`; a `col`, `lty` or `pch` given there wins and is
  recycled over the settings.

## Value

A list with `col`, `lty` and `pch`, each of length `k`.

## Details

A continuous family is separated by color and line type, a discrete one
by color and symbol. Dashing a stem is what a line type would do there,
and a dashed stem reads as a broken one: at a support of any size the
panel fills with fragments that cross each other. The symbol carries the
same information without drawing anything extra, the point at the top of
the stem being already there.
