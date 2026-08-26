# The Legend Entries and the Title of a Distribution Plot

Splits a set of plot settings into what varies and what does not: the
varying parameters become one legend entry per curve, and the fixed ones
are stated once in the title beside the family's name. A reader then
sees which parameter the panel is about without counting curves, and
reads the held values without a second key.

Where nothing varies there is one curve and no legend, and the title
carries every parameter. Where nothing is fixed the title is the
family's name alone. Values are rounded to three decimals.

## Usage

``` r
plot_labels(x, ps)
```

## Arguments

- x:

  A distribution object, read for `distrib_name` and `params`.

- ps:

  The value of
  [`plot_settings()`](https://statmodels7.github.io/distributions7/reference/plot_settings.md),
  read for `settings` and `varying`.

## Value

A list with two components: `legend`, a character vector of length
`ps$k` or `NULL` when nothing varies, and `main`, a character string of
length 1 carrying a newline between the family's name and the fixed
values.

## See also

[`plot_settings()`](https://statmodels7.github.io/distributions7/reference/plot_settings.md),
which supplies `ps`;
[`plot_legend_side()`](https://statmodels7.github.io/distributions7/reference/plot_legend_side.md)
for where the legend goes;
[`plot.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.continuous_distrib.md)
and
[`plot.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.discrete_distrib.md),
the callers.
