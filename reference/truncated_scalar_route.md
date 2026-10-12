# The Scalar Route of a Truncated Family

Returns the route of a
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
family, `"TruncatedDiscreteDistrib|<inner name>"` or
`"TruncatedContinuousDistrib|<inner name>"`, with the inner family's
constants followed by the truncation points `lower` and `upper`, as
[`distrib_scalar_route()`](https://statmodels7.github.io/distributions7/reference/distrib_scalar_route.md)
describes, or `NULL` when the parent has no route. A continuous parent
needs a family with a compiled distribution function and at most four
parameters, which may be wrapped by any chain of
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md),
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md)
and
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md)
(the last with one of the ready-made transformers).

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
