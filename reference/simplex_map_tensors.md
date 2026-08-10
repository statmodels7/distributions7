# Derivative Tensors of a Simplex-Valued Map

Collects \\\partial^{S}\mu_j\\ for every multiset \\S\\ of free indices
up to the requested order, keyed by the sorted tuple, as a list of
numeric vectors over the coordinates \\j\\.

## Usage

``` r
simplex_map_tensors(s, eta, order)
```

## Arguments

- s:

  A parameters7 `simplex` parameter.

- eta:

  The free vector.

- order:

  The highest order required, 1 to 4.

## Value

A named list of numeric vectors, keyed as `"1"`, `"1,2"` and so on.
