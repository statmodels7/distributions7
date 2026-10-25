# The Anchors of the Tail Correction of the Skew t

Returns the points at which
[`loc_scale_tail_weights()`](https://statmodels7.github.io/distributions7/reference/loc_scale_tail_weights.md)
places the weight of the tail past the last node of
[`loc_scale_rule()`](https://statmodels7.github.io/distributions7/reference/loc_scale_rule.md),
and the grid of the rule. There are nine anchors in each tail, equally
spaced in \\\log\|z\|\\ from \\10^{13}\\ to \\10^{150}\\.

## Usage

``` r
loc_scale_tail_anchors()
```

## Value

A list with the anchors `z` (positive; the left tail uses their
negatives), the step `h` of the rule, the constant `c0` of its
substitution and the index `kmax` of its last node.
