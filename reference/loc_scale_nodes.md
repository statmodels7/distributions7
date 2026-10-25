# The Nodes of the Location-Scale Rule for a Family

Returns the nodes and weights of
[`loc_scale_rule()`](https://statmodels7.github.io/distributions7/reference/loc_scale_rule.md),
extended for a skew t by the eighteen anchor nodes of its tail
correction, nine in each tail, with weight zero.
[`loc_scale_weights()`](https://statmodels7.github.io/distributions7/reference/loc_scale_weights.md)
gives the weights of the anchors for each degrees of freedom. For every
other family the rule is returned unchanged.

## Usage

``` r
loc_scale_nodes(distrib)
```

## Arguments

- distrib:

  A distribution object.

## Value

A list with the nodes `x` and the weights `w`.
