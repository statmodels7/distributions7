# Colors, Line Types and Symbols for a Set of Curves

Returns the three visual keys by which a graph distinguishes its
settings, and each key is recycled to length \\k\\. A single curve is
drawn in the main colour of the toolkit's graphs, solid and with a
filled circle. Several take a qualitative palette, the six line types
and six distinguishable symbols, so that a printed copy with no color is
still readable.

A continuous family draws with the color and the line type. A discrete
family draws with the color and the symbol, because a line type would
dash the stems and a dashed stem reads as a broken one.

## Usage

``` r
plot_keys(k)
```

## Arguments

- k:

  How many settings are drawn, a single positive integer.

## Value

A list with components `col`, `lty` and `pch`, each of length `k`.

## See also

[`plot.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.continuous_distrib.md)
and
[`plot.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.discrete_distrib.md),
which read different pairs of these;
[`plot_settings()`](https://statmodels7.github.io/distributions7/reference/plot_settings.md)
for `k`.
