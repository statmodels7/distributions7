# The Location and Scale Parts of the Response Curvature's Derivative

Computes \\\partial\ell^{(yy)}/\partial\mu = -B'/\sigma^3\\ and
\\\partial\ell^{(yy)}/\partial\sigma = -(zB' + 2B)/\sigma^3\\, with \\B
= \sigma^2\ell^{(yy)}\\ and \\B' = \sigma^3\ell^{(yyy)}\\, for a family
whose response enters only through \\z = (y-\mu)/\sigma\\. It is the
arithmetic behind every location-scale method of
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md).

## Usage

``` r
loc_scale_cross2_block(distrib, y, theta)
```

## Arguments

- distrib:

  An object inheriting from class `distrib`, whose first two parameters
  are a location and a scale in that order.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, location first and scale second.

## Value

A list of two numeric vectors, UNNAMED, in location-then-scale order.
The callers name them.

## Details

The same identity
[`loc_scale_cross_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross_block.md)
uses, one derivative further in the response: it reads the family's own
[`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
instead of differencing its
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md).
A shape parameter beyond the two is not covered, and
[`partial_loc_scale_cross2_y()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross2_y.md)
differences those components instead.

## Notation

\\\ell\\ is the log-density of one observation, \\\ell^{(y)}\\ and
\\\ell^{(yy)}\\ its first and second derivatives in the response, \\z =
(y-\mu)/\sigma\\ the standardized residual, and \\A\\ and \\B\\ the two
standardized quantities \\\sigma\ell^{(y)}\\ and
\\\sigma^2\ell^{(yy)}\\, whose derivatives are taken in \\z\\.

## See also

[`loc_scale_cross2_y()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross2_y.md)
and
[`partial_loc_scale_cross2_y()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross2_y.md),
the two method bodies built on it,
[`loc_scale_theta2_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_theta2_block.md)
for the next order in \\\theta\\, and
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
for the generic.

## Examples

``` r
d <- logistic_distrib()
y <- c(-0.7, 0.3, 1.4)
theta <- list(mu = 0.3, sigma = 1.2)

blk <- distributions7:::loc_scale_cross2_block(d, y, theta)
vapply(blk, function(z) z[1], numeric(1))
#> [1] 0.09632533 0.40854287

# The identity written out.
s <- 1.2
z <- (y - 0.3) / s
B <- s^2 * distrib_hess_y(d, y, theta)
B1 <- s^3 * distrib_deriv3_y(d, y, theta)
c(mu = (-B1 / s^3)[1], sigma = (-(z * B1 + 2 * B) / s^3)[1])
#>         mu      sigma 
#> 0.09632533 0.40854287 

# It is what the family's distrib_cross2_y method returns, named.
vapply(distrib_cross2_y(d, y, theta), function(z) z[1], numeric(1))
#>         mu      sigma 
#> 0.09632533 0.40854287 
```
