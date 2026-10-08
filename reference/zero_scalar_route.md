# The Scalar Route of a Zero Wrapper

Returns the route of a family built by
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
or
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md):
the name `"<class>|<inner name>"` and the inner family's constants, as
[`distrib_scalar_route()`](https://statmodels7.github.io/distributions7/reference/distrib_scalar_route.md)
describes, or `NULL` when the parent has no route of its own.

## Usage

``` r
zero_scalar_route(distrib, ...)
```

## Arguments

- distrib:

  A zero-inflated or zero-adjusted family.

- ...:

  Unused.

## Value

A list with components `name` and `constants`, or `NULL`.
