# Default Hyperparameter Hessian of the Response Gradient

Falls back to one tensor stencil on the highest quantity the family
implements itself, through
[`mixed_tensor_derivatives()`](https://statmodels7.github.io/distributions7/reference/mixed_tensor_derivatives.md):
a central difference of an analytic
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
in each parameter, a difference in the response of an analytic Hessian,
or one mixed stencil of the log-density, whichever leaves the smallest
order. Registering the fallback on `continuous_distrib` gives the
third-order mixed derivative to every continuous family, whether or not
it writes one out.

## Arguments

- distrib:

  A `continuous_distrib` object with no closed form of its own.

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

A named list with one numeric vector per unordered pair of parameters,
keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## Details

The stencil is never applied to a quantity that is itself a fallback.
Before 0.76.0 the fallback differenced
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
in the parameters whether or not that was analytic, and where it was a
fallback that differenced a difference in the same parameter.

## See also

[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
for the generic,
[`distrib_hess_y_hess.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.continuous_distrib.md)
for the fourth-order twin, and
[`numerical_theta2_y()`](https://statmodels7.github.io/distributions7/reference/numerical_theta2_y.md),
which does the work.

## Examples

``` r
# The gamma writes no closed form, so this method answers for it.
d <- gamma2_distrib()
y <- c(0.5, 1, 2)
theta <- list(mu = 2, sigma2 = 1)
attr(S7::method(distrib_grad_y_hess, S7::S7_class(d)), "signature")[[1]]
#> <distributions7::continuous_distrib> class
#> @ parent     : <distributions7::distrib>
#> @ constructor: function(distrib_name, dimension, bounds, params, params_interpretation, n_params, params_bounds, link_params, params_smooth) {...}
#> @ validator  : <NULL>
#> @ properties :
#>  $ distrib_name         : <character>          
#>  $ dimension            : <character>          
#>  $ bounds               : <integer> or <double>
#>  $ params               : <character>          
#>  $ params_interpretation: <character>          
#>  $ n_params             : <integer> or <double>
#>  $ params_bounds        : <list>               
#>  $ link_params          : <list>               
#>  $ params_smooth        : <logical>            

g <- distrib_grad_y_hess(d, y, theta)
c(g$mu_mu[1], g$mu_sigma2[1], g$sigma2_sigma2[1])
#> [1]  4 -7 12

# Against a numerical Hessian of the response gradient.
f <- function(v) distrib_grad_y(d, y[1], list(mu = v[1], sigma2 = v[2]))
numDeriv::hessian(f, c(2, 1))
#>      [,1] [,2]
#> [1,]    4   -7
#> [2,]   -7   12
```
