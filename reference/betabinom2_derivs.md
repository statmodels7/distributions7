# Every Component of a Beta-Binomial Derivative in the Shapes

Assembles the full set of derivative components of a given order by
calling
[`betabinom2_component()`](https://statmodels7.github.io/distributions7/reference/betabinom2_component.md)
once per distinct multi-index, and names them as
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
names them.

## Usage

``` r
betabinom2_derivs(y, a, b, n, order, params)
```

## Arguments

- y:

  A numeric vector of counts.

- a, b:

  The two shapes, each a numeric vector of length 1 or of the length of
  `y`, strictly positive.

- n:

  The size, a single positive integer.

- order:

  The derivative order, an integer from 1 to 4.

- params:

  The parameter names, `c("alpha", "beta")` for this family, used to
  build the component names and the multi-indices in the same order.

## Value

A named list of component vectors, one per distinct multi-index of the
given order: two at order 1, three at order 2, four at order 3 and five
at order 4. Each has the recycled length of the inputs.

## See also

[`betabinom2_component()`](https://statmodels7.github.io/distributions7/reference/betabinom2_component.md)
for one component,
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
and
[`deriv_indices()`](https://statmodels7.github.io/distributions7/reference/deriv_indices.md)
for the enumeration, and
[`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md)
for the family.
