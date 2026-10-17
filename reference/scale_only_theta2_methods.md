# The Mixed Grid of a Family With No Location

Builds the
[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
and
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md)
method bodies for a family with a scale and no location: the scale's own
pair comes from
[`scale_only_theta2()`](https://statmodels7.github.io/distributions7/reference/scale_only_theta2.md)
in closed form, and every pair that involves a shape parameter comes
from one central difference of the first-order component through
[`numerical_theta2_y()`](https://statmodels7.github.io/distributions7/reference/numerical_theta2_y.md).
In the fourth-order body that component, from
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md),
is itself differenced for the shape parameter.

## Usage

``` r
scale_only_theta2_methods(at = 1L)
```

## Arguments

- at:

  The index of the scale parameter in `distrib@params`. The default is
  `1`.

## Value

A list of two method bodies, named `grad2` and `hess2`, each with the
generics' signature.

## Details

A one-parameter family returns only the closed component and differences
nothing. The bodies are registered on the exponential, which is that
case, and on the Weibull, whose second parameter is a shape. The
generalized Pareto has compiled methods of its own.

## Notation

\\\ell\\ is the log-density of one observation, \\\ell^{(y)}\\ and
\\\ell^{(yy)}\\ its first and second derivatives in the response,
\\\sigma\\ the scale and \\z\\ the standardized response, \\y/\sigma\\
where there is no location and \\(y-\mu)/\sigma\\ where there is.

## See also

[`scale_only_theta2()`](https://statmodels7.github.io/distributions7/reference/scale_only_theta2.md)
for the closed component,
[`numerical_theta2_y()`](https://statmodels7.github.io/distributions7/reference/numerical_theta2_y.md)
for the differenced ones, and
[`partial_theta2()`](https://statmodels7.github.io/distributions7/reference/partial_theta2.md),
the location-scale counterpart.

## Examples

``` r
y <- c(0.4, 1.1, 2.3)

# For one parameter, the method returns only the closed component.
d <- exponential_distrib()
distrib_grad_y_hess(d, y, list(mu = 1.5))
#> $mu_mu
#> [1] -0.5925926 -0.5925926 -0.5925926
#> 

# For two parameters, the scale pair is closed and the two pairs that
# involve the shape are differenced.
dw <- weibull1_distrib()
theta <- list(mu = 1.5, sigma = 1.3)
g <- distrib_grad_y_hess(dw, y, theta)
names(g)
#> [1] "mu_mu"       "sigma_sigma" "mu_sigma"   
vapply(g, function(z) z[1], numeric(1))
#>       mu_mu sigma_sigma    mu_sigma 
#>  -0.7746961   0.1669799   0.1094877 

# The result is compared with a numerical Hessian of the response gradient.
f <- function(v) distrib_grad_y(dw, y[1], list(mu = v[1], sigma = v[2]))
numDeriv::hessian(f, c(1.5, 1.3))
#>            [,1]      [,2]
#> [1,] -0.7746961 0.1094877
#> [2,]  0.1094877 0.1669799

# The scale pair is exactly the closed component.
identical(g$mu_mu,
          distributions7:::scale_only_theta2(dw, y, theta, 1L, 1L))
#> [1] TRUE
```
