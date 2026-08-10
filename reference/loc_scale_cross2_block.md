# The Location and Scale Components of the Response Curvature's Derivative

\\\partial\ell^{(yy)}/\partial\mu = -B'/\sigma^3\\ and
\\\partial\ell^{(yy)}/\partial\sigma = -(zB' + 2B)/\sigma^3\\, with \\B
= \sigma^2\ell^{(yy)}\\ and \\B' = \sigma^3\ell^{(yyy)}\\, for a family
whose response enters only through \\z = (y-\mu)/\sigma\\.

## Usage

``` r
loc_scale_cross2_block(distrib, y, theta)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, location first and scale second.

## Value

A list of two component vectors, unnamed.

## Details

The same identity as
[`loc_scale_cross_block`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross_block.md)
one derivative further in the response, and it reads the family's own
third response derivative rather than differencing its second.

## See also

[`distrib_cross2_y`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
