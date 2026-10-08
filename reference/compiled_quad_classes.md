# Families With a Compiled CDF Quadrature

`compiled_quad_classes` names the classes whose compiled distribution
function integrates its shape derivatives, which
[`compiled_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf_deriv_k.md)
serves at orders three and four; `compiled_ls_classes` names those among
them that are location-scale in their first two parameters.

## Usage

``` r
compiled_quad_classes

compiled_ls_classes
```

## Format

An object of class `character` of length 14.

An object of class `character` of length 5.

## Value

Character vectors of S7 class names.
