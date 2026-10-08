# The Scalar Route of a Truncated Family

Returns the route of a
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
family, `"TruncatedDiscreteDistrib|<inner name>"` or
`"TruncatedContinuousDistrib|<inner name>"`, with the inner family's
constants followed by the truncation points `lower` and `upper`, as
[`distrib_scalar_route()`](https://statmodels7.github.io/distributions7/reference/distrib_scalar_route.md)
describes, or `NULL` when the parent has no route. A continuous parent
needs a compiled distribution function, and the only wrapper it may
carry is
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md).

## Usage

``` r
truncated_scalar_route(distrib, ...)
```

## Arguments

- distrib:

  A truncated family.

- ...:

  Unused.

## Value

A list with components `name` and `constants`, or `NULL`.
