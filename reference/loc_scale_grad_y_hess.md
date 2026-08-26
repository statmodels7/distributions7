# Second-Order Mixed Derivatives of a Location-Scale Family

The
[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
and
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md)
method bodies shared by the families that are location-scale in BOTH
their parameters. Each one calls
[`loc_scale_theta2_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_theta2_block.md)
at its own order and returns the three components in
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)'s
order, so neither derives anything of its own.

## Usage

``` r
loc_scale_grad_y_hess(distrib, y, theta, scale = c("parameter", "link"), ...)

loc_scale_hess_y_hess(distrib, y, theta, scale = c("parameter", "link"), ...)
```

## Arguments

- distrib:

  An object inheriting from class `distrib`, whose first two parameters
  are a location and a scale in that order.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- scale:

  One of `"parameter"` or `"link"`, applied by the generic before
  dispatch.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of three numeric vectors, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## Details

Registered on the logistic, Cauchy, Gumbel and Laplace. A family
carrying a shape parameter beyond the two takes
[`partial_loc_scale_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_grad_y_hess.md)
instead, where the pairs touching that shape are differenced.

Both bodies read the family's own third and fourth response derivatives,
so a family supplying
[`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
and
[`distrib_deriv4_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
gets these in closed form with no further algebra.

## Notation

\\\ell\\ is the log-density of one observation, \\\ell^{(y)}\\ and
\\\ell^{(yy)}\\ its first and second derivatives in the response, \\z =
(y-\mu)/\sigma\\ the standardized residual, and \\A\\ and \\B\\ the two
standardized quantities \\\sigma\ell^{(y)}\\ and
\\\sigma^2\ell^{(yy)}\\, whose derivatives are taken in \\z\\.

## See also

[`loc_scale_theta2_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_theta2_block.md),
which does the work,
[`partial_loc_scale_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_grad_y_hess.md)
for a family with a shape parameter, and
[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
and
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md)
for the generics.

## Examples

``` r
d <- gumbel_distrib()
y <- c(-0.7, 0.3, 1.4)
theta <- list(mu = 0.3, sigma = 1.2)

g3 <- distrib_grad_y_hess(d, y, theta)
vapply(g3, function(z) z[1], numeric(1))
#>       mu_mu sigma_sigma    mu_sigma 
#>    1.331583    6.869081   -3.772819 

# Against a numerical Hessian of the response gradient.
f <- function(v) distrib_grad_y(d, y[1], list(mu = v[1], sigma = v[2]))
numDeriv::hessian(f, c(0.3, 1.2))
#>           [,1]      [,2]
#> [1,]  1.331583 -3.772819
#> [2,] -3.772819  6.869081

# The fourth order takes the same body one order up.
g4 <- distrib_hess_y_hess(d, y, theta)
vapply(g4, function(z) z[1], numeric(1))
#>       mu_mu sigma_sigma    mu_sigma 
#>   -1.109653  -12.976772    4.253669 
h <- function(v) distrib_hess_y(d, y[1], list(mu = v[1], sigma = v[2]))
numDeriv::hessian(h, c(0.3, 1.2))
#>           [,1]       [,2]
#> [1,] -1.109653   4.253669
#> [2,]  4.253669 -12.976772
```
