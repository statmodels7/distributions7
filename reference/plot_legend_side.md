# Where to Put the Key

The corner of the panel a legend goes in: the one the curves leave
emptier.

## Usage

``` r
plot_legend_side(y, dens)
```

## Arguments

- y:

  The evaluation points.

- dens:

  A list of densities, one per setting.

## Value

`"topright"` or `"topleft"`.

## Details

A right-skewed family puts its mass on the left, so a key fixed at the
top right never meets it, while a left-skewed one is covered by exactly
that choice. The side is therefore read off the drawing: the density's
center of mass is compared with the middle of the horizontal range, and
the key goes to the other side.
