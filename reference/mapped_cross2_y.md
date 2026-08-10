# Second-Response Mixed Derivatives Through a Map

The same first-order chain rule
[`mapped_cross_y`](https://statmodels7.github.io/distributions7/reference/mapped_cross_y.md)
takes, on the second response derivative.

## Usage

``` r
mapped_cross2_y(distrib, parent, th_par, maps, y)
```

## Arguments

- distrib:

  The distribution in the new parametrization.

- parent:

  The parent distribution.

- th_par:

  The parent's parameters at the new ones.

- maps:

  The map's keyed partial tables.

- y:

  A numeric vector of observations.

## Value

A named list with one numeric vector per new parameter.

## See also

[`mapped_cross_y`](https://statmodels7.github.io/distributions7/reference/mapped_cross_y.md)
