# Does a Reparametrization Carry the Map's Closed Derivatives?

Returns `TRUE` when the map's partial derivatives of a reparametrized
family were supplied through `map_derivs`, and `FALSE` when they come
from the finite-difference tables of
[`reparam_stencil_derivs()`](https://statmodels7.github.io/distributions7/reference/reparam_stencil_derivs.md).

## Usage

``` r
reparam_has_map_derivs(distrib)
```

## Arguments

- distrib:

  A reparametrized family.

## Value

A single logical.

## Details

The numerical tables lose accuracy quickly with the order, so a chain of
the third or fourth order that reads them is no better than the stencil
it would replace. The two routes are told apart by the body of the
`reparam_derivs` function, which is the one of
[`reparam_stencil_derivs()`](https://statmodels7.github.io/distributions7/reference/reparam_stencil_derivs.md)
in the numerical case.
