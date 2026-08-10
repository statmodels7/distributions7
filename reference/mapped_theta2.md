# A Second-Order Mixed Derivative Through a Map

The parent's paired components carried onto a new parametrization by the
second-order chain rule.

## Usage

``` r
mapped_theta2(distrib, parent, th_par, maps, y, first, second)
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

- first:

  The parent's first-order components, keyed by parameter.

- second:

  The parent's second-order components, keyed by pair.

## Value

A named list keyed by the new parameters' pairs.

## See also

[`mapped_cross_y`](https://statmodels7.github.io/distributions7/reference/mapped_cross_y.md),
[`reparam_tables`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md)
