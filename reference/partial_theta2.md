# Splice the Closed Location-Scale Pairs Into the Fallback

At order 1, differences every component of
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
through
[`numerical_theta2_y()`](https://statmodels7.github.io/distributions7/reference/numerical_theta2_y.md)
and then replaces the three that
[`loc_scale_theta2_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_theta2_block.md)
gives in closed form. At order 2, a pair of a location or scale
parameter with a shape is one central difference, in the shape, of the
closed component of
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md),
and a pair of shapes is one stencil on
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md),
a second-order stencil for a repeated shape and a mixed one for two
shapes; the shape components of
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
are not differenced again. It is what
[`partial_loc_scale_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_grad_y_hess.md)
and
[`partial_loc_scale_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_grad_y_hess.md)
both run, at order 1 and 2 respectively.

## Usage

``` r
partial_theta2(distrib, y, theta, order)
```

## Arguments

- distrib:

  An object inheriting from class `distrib`, whose first two parameters
  are a location and a scale.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- order:

  `1` to splice into the derivative of
  [`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
  `2` into that of
  [`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md).
  It selects both the first-order reference differenced,
  [`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
  or
  [`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md),
  and the order of the closed block.

## Value

A named list with one numeric vector per unordered pair of parameters,
keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## Details

[`numerical_theta2_y()`](https://statmodels7.github.io/distributions7/reference/numerical_theta2_y.md)
takes one central difference per parameter, and each difference yields
every component that involves that parameter, so replacing the three
location-scale pairs afterwards adds no evaluation of the first-order
function. The result is keyed by
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
and the closed block's three names index into it directly.

## Notation

\\\ell\\ is the log-density of one observation, \\\ell^{(y)}\\ and
\\\ell^{(yy)}\\ its first and second derivatives in the response, \\z =
(y-\mu)/\sigma\\ the standardized residual, and \\A\\ and \\B\\ the two
standardized quantities \\\sigma\ell^{(y)}\\ and
\\\sigma^2\ell^{(yy)}\\, whose derivatives are taken in \\z\\.

## See also

[`loc_scale_theta2_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_theta2_block.md)
for the closed three,
[`numerical_theta2_y()`](https://statmodels7.github.io/distributions7/reference/numerical_theta2_y.md)
for the differenced rest, and
[`partial_loc_scale_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_grad_y_hess.md),
its caller.

## Examples

``` r
d <- skewt_distrib()
y <- c(-0.7, 0.3, 1.4)
theta <- list(mu = 0.3, sigma = 1.2, alpha = 0.7, nu = 6)

out <- distributions7:::partial_theta2(d, y, theta, 1L)
names(out)
#>  [1] "mu_mu"       "sigma_sigma" "alpha_alpha" "nu_nu"       "mu_sigma"   
#>  [6] "mu_alpha"    "mu_nu"       "sigma_alpha" "sigma_nu"    "alpha_nu"   

# The three location-scale pairs are spliced in unchanged.
cl <- distributions7:::loc_scale_theta2_block(d, y, theta, 1L)
all(vapply(names(cl), function(k) identical(out[[k]], cl[[k]]), logical(1)))
#> [1] TRUE

# And it is exactly what the family's method returns.
identical(out, distrib_grad_y_hess(d, y, theta))
#> [1] TRUE
```
