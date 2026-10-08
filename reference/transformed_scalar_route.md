# The Scalar Route of a Transformed Family

Returns the route of a family built by
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md):
the name `"TransformedDistrib:<code>|<inner name>"` and the inner
family's constants followed by the transformer's parameters, as
[`distrib_scalar_route()`](https://statmodels7.github.io/distributions7/reference/distrib_scalar_route.md)
describes, or `NULL` when the parent has no route of its own or the
transformer is not one of the ready-made ones.

## Usage

``` r
transformed_scalar_route(distrib, ...)
```

## Arguments

- distrib:

  A transformed family.

- ...:

  Unused.

## Value

A list with components `name` and `constants`, or `NULL`.
