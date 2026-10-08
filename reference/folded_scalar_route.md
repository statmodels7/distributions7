# The Scalar Route of a Folded Family

Returns the route of a
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md)
family, `"FoldedDistrib|<inner name>"` with the inner family's
constants, or `NULL` when the parent has no route or is not a
location-scale family on the real line (possibly with some parameters
fixed by
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md)),
the expected information's quadrature reading the parent's center and
scale.

## Usage

``` r
folded_scalar_route(distrib, ...)
```

## Arguments

- distrib:

  A folded family.

- ...:

  Unused.

## Value

A list with components `name` and `constants`, or `NULL`.
