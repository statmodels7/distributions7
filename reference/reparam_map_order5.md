# The Map's Partial Derivatives of the Fifth Order

Adds to the keyed tables of a reparametrization the partials of the
fifth order, each the derivative along its largest coordinate of the
partial of the fourth order below it.

## Usage

``` r
reparam_map_order5(distrib, theta, maps)
```

## Arguments

- distrib:

  A reparametrized family whose `reparam_derivs` supplies the fourth
  order.

- theta:

  A named list of the new parameters.

- maps:

  The tables to the fourth order, from
  [`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md).

## Value

`maps` with the partials of the fifth order added under their keys.

## Details

The contract of `map_derivs` stops at the fourth order, and the fourth
derivative of the expected information reads the fifth. A partial of the
fifth order is the derivative, along one new coordinate \\\phi_c\\, of
an analytic partial of the fourth order. It is formed by a central
difference of the tables at \\\phi_c \pm h\\ and at \\\pm h/2\\,
extrapolated as \\(4D(h/2) - D(h))/3\\, with \\h\\ the relative step
\\2\cdot 10^{-4}\\ of the coordinate kept inside its bounds by
[`fd_steps()`](https://statmodels7.github.io/distributions7/reference/fd_steps.md).
The tables are closed forms evaluated at four shifted points per
coordinate, with no evaluation of the parent, and the difference is
taken once on an analytic quantity.
