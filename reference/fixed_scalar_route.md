# The Scalar Route of a Fixed Family

Returns the route of a family built by
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md):
the name `"<class>:<mask>|<inner name>"` and the inner family's
constants followed by the fixed values, as
[`distrib_scalar_route()`](https://statmodels7.github.io/distributions7/reference/distrib_scalar_route.md)
describes, or `NULL` when the parent has no route of its own or a fixed
value varies by observation.

## Usage

``` r
fixed_scalar_route(distrib, ...)
```

## Arguments

- distrib:

  A fixed family.

- ...:

  Unused.

## Value

A list with components `name` and `constants`, or `NULL`.
