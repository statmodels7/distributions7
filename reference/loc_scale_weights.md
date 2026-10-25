# The Weights of the Location-Scale Rule for Each Shape

Returns the weights of the nodes of
[`loc_scale_nodes()`](https://statmodels7.github.io/distributions7/reference/loc_scale_nodes.md)
for each row of the shapes `U`, as one vector of length
`length(rule$x) * nrow(U)`. For a skew t the weights of the anchor nodes
are those of
[`loc_scale_tail_weights()`](https://statmodels7.github.io/distributions7/reference/loc_scale_tail_weights.md),
so that the integral past the last node of the rule is included; for
every other family they are the weights of the rule.

## Usage

``` r
loc_scale_weights(distrib, U, rule = loc_scale_nodes(distrib))
```

## Arguments

- distrib:

  A distribution object.

- U:

  The matrix of distinct shapes, one row per shape.

- rule:

  The list returned by
  [`loc_scale_nodes()`](https://statmodels7.github.io/distributions7/reference/loc_scale_nodes.md).

## Value

A numeric vector.
