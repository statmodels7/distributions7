# Labels for a Plot of Several Settings

The legend entries and the title of a plot drawing several settings: the
parameters that vary go in the legend, one entry per curve, and those
held fixed go in the title, where they are stated once.

## Usage

``` r
plot_labels(x, ps)
```

## Arguments

- x:

  A distribution object.

- ps:

  The value of
  [`plot_settings`](https://statmodels7.github.io/distributions7/reference/plot_settings.md).

## Value

A list with `legend` (character, length `ps$k`, or `NULL` when nothing
varies) and `main`.
