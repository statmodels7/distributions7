# The Mixed Grid of a Family Written Out Against a Tabulated Map

`distrib_cross2_y`, `distrib_grad_y_hess` and `distrib_hess_y_hess` for
a family that carries its own class and kernels while its map onto a
parent is tabulated.

## Usage

``` r
mapped_theta2_methods(parent, th_par, tables)
```

## Arguments

- parent:

  The parent distribution.

- th_par:

  A function of the new parameters returning the parent's.

- tables:

  A function of the new parameters returning the map's keyed partial
  tables.

## Value

A list of three method bodies, named `cross2`, `grad2` and `hess2`.

## See also

[`mapped_theta2`](https://statmodels7.github.io/distributions7/reference/mapped_theta2.md),
[`reparam_map_derivs`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md)
