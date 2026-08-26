# Colors, Line Types and Symbols for a Set of Curves

Returns the three visual keys a plot method distinguishes its settings
by, each recycled to length \\k\\. A single curve is black, solid and a
filled circle. Several take a qualitative palette, the six base line
types and six distinguishable symbols, so that a printed copy with no
color is still readable.

A continuous family draws with the color and the line type; a discrete
one draws with the color and the **symbol**, because dashing a stem is
what a line type does there and a dashed stem reads as a broken one.

## Usage

``` r
plot_keys(k, dots = list())
```

## Arguments

- k:

  How many settings are drawn, a single positive integer.

- dots:

  The caller's `...`, as a list. A `col`, `lty` or `pch` given there
  wins and is recycled over the settings; anything else is ignored.

## Value

A list with components `col`, `lty` and `pch`, each of length `k`.

## See also

[`plot.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.continuous_distrib.md)
and
[`plot.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.discrete_distrib.md),
which read different pairs of these;
[`plot_settings()`](https://statmodels7.github.io/distributions7/reference/plot_settings.md)
for `k`.
