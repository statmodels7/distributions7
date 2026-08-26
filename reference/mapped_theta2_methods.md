# The Mixed Grid of a Family Written Out Against a Tabulated Map

Builds the
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md),
[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
and
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md)
method bodies for a family that carries its own class and kernels while
its map onto a parent is tabulated. Each body is the chain rule of
[`mapped_cross2_y()`](https://statmodels7.github.io/distributions7/reference/mapped_cross2_y.md)
or
[`mapped_theta2()`](https://statmodels7.github.io/distributions7/reference/mapped_theta2.md)
with the parent, the map and the partial tables closed over, so no new
algebra is involved.

## Usage

``` r
mapped_theta2_methods(parent, th_par, tables)
```

## Arguments

- parent:

  The parent distribution, built once and closed over.

- th_par:

  A function of the new parameters returning the parent's, as a named
  list.

- tables:

  A function of the new parameters returning the map's keyed partial
  tables, one of the `md_*` functions of `reparam_maps.R`.

## Value

A list of three method bodies, named `cross2`, `grad2` and `hess2`, each
with the generics' signature.

## Details

The families it serves have their own classes and their own compiled
kernels, so they are NOT `ReparamContinuousDistrib` and do not inherit
its methods; what they share with it is the map, which `reparam_maps.R`
already tabulates. Registered on `gaussian2`, `gaussian3`, `laplace2`
and `invgauss2`.

## Notation

\\\ell\\ is the log-density of one observation, \\\ell^{(y)}\\ and
\\\ell^{(yy)}\\ its first and second derivatives in the response,
\\\sigma\\ the scale and \\z\\ the standardized response, \\y/\sigma\\
where there is no location and \\(y-\mu)/\sigma\\ where there is.

## See also

[`mapped_theta2()`](https://statmodels7.github.io/distributions7/reference/mapped_theta2.md)
and
[`mapped_cross2_y()`](https://statmodels7.github.io/distributions7/reference/mapped_cross2_y.md),
which do the work,
[`reparam_map_derivs()`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md)
for the tables, and
[`scale_only_theta2_methods()`](https://statmodels7.github.io/distributions7/reference/scale_only_theta2_methods.md)
for the other shape this file covers.

## Examples

``` r
# gaussian2 carries (mu, sigma2) against a parent in (mu, sigma).
d <- gaussian2_distrib()
y <- c(-0.7, 0.3, 1.4)
theta <- list(mu = 0.3, sigma2 = 1.44)

g3 <- distrib_grad_y_hess(d, y, theta)
vapply(g3, function(z) z[1], numeric(1))
#>         mu_mu sigma2_sigma2     mu_sigma2 
#>     0.0000000     0.6697960    -0.4822531 

# Against a numerical Hessian of the response gradient in the family's own
# parameters, which shares none of the map's arithmetic.
f <- function(v) distrib_grad_y(d, y[1], list(mu = v[1], sigma2 = v[2]))
numDeriv::hessian(f, c(0.3, 1.44))
#>               [,1]       [,2]
#> [1,]  9.210739e-13 -0.4822531
#> [2,] -4.822531e-01  0.6697960

# laplace2 carries the rate, so its map is 1 / lambda.
dl <- laplace2_distrib()
vapply(distrib_cross2_y(dl, c(-0.7, 1.1, 1.4), list(mu = 0.3, lambda = 0.8)),
       function(z) z[1], numeric(1))
#>     mu lambda 
#>      0      0 
```
