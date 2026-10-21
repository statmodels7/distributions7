# Placing the Legend

Returns the corner of the panel in which the legend is drawn, the one
that the curves leave emptier. A right-skewed family puts its mass on
the left, so a key at the top right does not meet it; a left-skewed
family puts its mass on the right, where a key at the top right would
cover it.

The side is read off the drawing. The density's center of mass across
all the settings is compared with the middle of the horizontal range,
and the key goes to the other side.

## Usage

``` r
plot_legend_side(y, dens)
```

## Arguments

- y:

  The evaluation points, a numeric vector.

- dens:

  A list of densities, one per setting, in which each has the length of
  `y`. Non-finite and negative entries are treated as zero.

## Value

`"topright"` or `"topleft"`, a character string of length 1. The
function returns `"topright"` when no setting has any positive mass.

## See also

[`plot.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.continuous_distrib.md)
and
[`plot.discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.discrete_distrib.md),
the callers;
[`plot_labels()`](https://statmodels7.github.io/distributions7/reference/plot_labels.md)
for what goes in the key.
