# Mixed Derivatives Through a Map

Carries the parent's mixed block onto a new parametrization by the
first-order chain rule. That rule is all that is needed, because the
response derivative does not interact with a reparametrization of the
parameters.

## Usage

``` r
mapped_cross_y(distrib, parent, th_par, maps, y)
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

[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md),
[`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md)
