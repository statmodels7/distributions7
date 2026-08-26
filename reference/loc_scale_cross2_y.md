# Second-Response Mixed Derivatives of a Location-Scale Family

The
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
method bodies shared by the families whose response enters only through
\\z = (y-\mu)/\sigma\\. `loc_scale_cross2_y()` serves a family that is
location-scale in BOTH its parameters, and
`partial_loc_scale_cross2_y()` one carrying shape parameters beyond the
two: the location and scale components come from
[`loc_scale_cross2_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross2_block.md)
in closed form, and each shape component from one central difference of
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
through
[`numerical_cross2_y()`](https://statmodels7.github.io/distributions7/reference/numerical_cross2_y.md).

## Usage

``` r
loc_scale_cross2_y(distrib, y, theta, scale = c("parameter", "link"), ...)

partial_loc_scale_cross2_y(
  distrib,
  y,
  theta,
  scale = c("parameter", "link"),
  ...
)
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

A named list with one numeric vector per parameter, keyed by
`distrib@params`.

## Details

Registered on the logistic, Cauchy, Gumbel and Laplace for the first
body, and on the pseudo-Huber, skew normal and skew t for the second.
The gaussian and the Student t write their own out, and every remaining
continuous family takes
[`distrib_cross2_y.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.continuous_distrib.md).

## Notation

\\\ell\\ is the log-density of one observation, \\\ell^{(y)}\\ and
\\\ell^{(yy)}\\ its first and second derivatives in the response, \\z =
(y-\mu)/\sigma\\ the standardized residual, and \\A\\ and \\B\\ the two
standardized quantities \\\sigma\ell^{(y)}\\ and
\\\sigma^2\ell^{(yy)}\\, whose derivatives are taken in \\z\\.

## See also

[`loc_scale_cross2_block()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cross2_block.md)
for the closed components,
[`numerical_cross2_y()`](https://statmodels7.github.io/distributions7/reference/numerical_cross2_y.md)
for the differenced ones, and
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
for the generic.

## Examples

``` r
y <- c(-0.7, 0.3, 1.4)

# Fully location-scale: both components closed.
d <- cauchy_distrib()
theta <- list(mu = 0.3, sigma = 1.2)
vapply(distrib_cross2_y(d, y, theta), function(z) z[1], numeric(1))
#>         mu      sigma 
#>  0.9141734 -0.5154616 

# Against a numerical derivative of the response Hessian.
f <- function(v) distrib_hess_y(d, y[1], list(mu = v[1], sigma = v[2]))
numDeriv::grad(f, c(0.3, 1.2))
#> [1]  0.9141734 -0.5154616

# With a shape parameter, the third component is differenced and the first
# two are not.
ds <- skewnormal1_distrib()
ts <- list(mu = 0.3, sigma = 1.2, alpha = 0.7)
vapply(distrib_cross2_y(ds, y, ts), function(z) z[1], numeric(1))
#>          mu       sigma       alpha 
#> -0.03060209  1.60525989 -0.76774711 
g <- function(v) distrib_hess_y(ds, y[1], list(mu = v[1], sigma = v[2],
                                               alpha = v[3]))
numDeriv::grad(g, c(0.3, 1.2, 0.7))
#> [1] -0.03060209  1.60525989 -0.76774711
```
