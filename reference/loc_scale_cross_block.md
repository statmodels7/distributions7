# The Location and Scale Components of a Mixed Derivative

\\\partial^2\ell/\partial y\partial\mu = -\ell^{(yy)}\\ and
\\\partial^2\ell/\partial y\partial\sigma = -z\ell^{(yy)} -
\ell^{(y)}/\sigma\\, for a family whose response enters only through \\z
= (y-\mu)/\sigma\\.

## Usage

``` r
loc_scale_cross_block(distrib, y, theta)
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

## See also

[`distrib_cross_y`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
