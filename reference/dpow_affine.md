# Derivatives of a Power of an Affine Argument

Returns \\d^{k}u^{p}/dv^{k}\\ for \\u = v\\ or \\u = 1 - v\\, the two
shapes the elastic net's scale and its argument are built from.

## Usage

``` r
dpow_affine(u, p, k, inner)
```

## Arguments

- u:

  The base, already evaluated.

- p:

  The exponent.

- k:

  The derivative order.

- inner:

  The derivative of the base in the variable, 1 or -1.

## Value

A numeric vector.
