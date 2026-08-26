# One Parameter Setting Per Curve

Turns the `theta` a plot method was given into the list of scalar
settings it draws, one per curve. A component of length \\k \> 1\\ asks
for \\k\\ curves, so `list(mu = 0, sigma = c(1, 2, 4))` becomes three
settings sharing a mean and differing in scale.

## Usage

``` r
plot_settings(x, theta)
```

## Arguments

- x:

  A distribution object, read for `params`.

- theta:

  A named list, or a named numeric vector, holding one entry per
  parameter. A missing parameter signals an error naming what was
  expected; so does a component of length zero.

## Value

A list with three components: `settings`, a list of \\k\\ named
parameter lists each holding scalars; `k`, the number of curves; and
`varying`, a character vector naming the parameters that differ between
settings, which may be empty.

## Details

Every component must have length 1 or the same \\k\\. A length that
merely **divides** \\k\\ signals an error: R would recycle it without a
word, and a partial setting is far likelier to be a mistake than a
request.

The meaning is available only because a plot has no data to recycle a
parameter against.
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
and every derivative generic read a vector component as one value per
observation, which is why each setting is handed to them as scalars.

## See also

[`plot.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.continuous_distrib.md)
and
[`plot.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.discrete_distrib.md),
the callers;
[`plot_labels()`](https://statmodels7.github.io/distributions7/reference/plot_labels.md),
which turns `varying` into a legend;
[`plot_keys()`](https://statmodels7.github.io/distributions7/reference/plot_keys.md)
for the colors and symbols.
