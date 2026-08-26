# The Scale Component of the Response Curvature's Derivative

Builds the
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
method body for a family with a scale and no location. The scale's
component is the identity
[`loc_scale_cross2_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross2_block.md)
uses, read at \\z = y/\sigma\\: \\-(zB' + 2B)/\sigma^3\\ with \\B =
\sigma^2\ell^{(yy)}\\ and \\B' = \sigma^3\ell^{(yyy)}\\. Any shape
component is differenced through
[`numerical_cross2_y()`](https://statmodels7.github.io/distributions7/reference/numerical_cross2_y.md).

## Usage

``` r
scale_only_cross2_method(at = 1L)
```

## Arguments

- at:

  The index of the scale parameter in `distrib@params`. Default `1`.

## Value

A method body with
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)'s
signature, returning a named list with one numeric vector per parameter.

## Details

A one-parameter family returns the closed component alone. Registered on
the exponential, the Weibull and the generalized Pareto.

## Notation

\\\ell\\ is the log-density of one observation, \\\ell^{(y)}\\ and
\\\ell^{(yy)}\\ its first and second derivatives in the response,
\\\sigma\\ the scale and \\z\\ the standardized response, \\y/\sigma\\
where there is no location and \\(y-\mu)/\sigma\\ where there is.

## See also

[`loc_scale_cross2_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross2_block.md),
the location-scale form,
[`scale_only_theta2_methods()`](https://statmodels7.github.io/distributions7/reference/scale_only_theta2_methods.md)
for the next order in \\\theta\\, and
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
for the generic.

## Examples

``` r
y <- c(0.4, 1.1, 2.3)

# The generalized Pareto: a scale and a shape.
d <- gpd_distrib()
theta <- list(sigma = 1.5, xi = 0.3)
vapply(distrib_cross2_y(d, y, theta), function(z) z[1], numeric(1))
#>      sigma         xi 
#> -0.1834635  0.5362778 

# Against a numerical derivative of the response Hessian.
f <- function(v) distrib_hess_y(d, y[1], list(sigma = v[1], xi = v[2]))
numDeriv::grad(f, c(1.5, 0.3))
#> [1] -0.1834635  0.5362778

# The exponential has one parameter, so nothing is differenced.
distrib_cross2_y(exponential_distrib(), y, list(mu = 1.5))
#> $mu
#> [1] 0 0 0
#> 
```
