# Default Mixed Second-Response Derivatives for Continuous Distributions

Falls back to one central difference of
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
in each parameter, through
[`numerical_cross2_y()`](https://statmodels7.github.io/distributions7/reference/numerical_cross2_y.md).
Registering it on `continuous_distrib` gives the quantity to every
continuous family, whether or not it writes one out.

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

A named list with one numeric vector per parameter, keyed by
`distrib@params`.

## See also

[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
for the generic,
[`numerical_cross2_y()`](https://statmodels7.github.io/distributions7/reference/numerical_cross2_y.md),
which does the work, and
[`distrib_cross2_y.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.Gaussian1Distrib.md)
for a family that overrides it.

## Examples

``` r
# The gamma writes no closed form, so this method answers for it.
d <- gamma2_distrib()
y <- c(0.5, 1, 2)
theta <- list(mu = 2, sigma2 = 1)
attr(S7::method(distrib_cross2_y, S7::S7_class(d)), "signature")[[1]]
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

vapply(distrib_cross2_y(d, y, theta), function(z) z[1], numeric(1))
#>     mu sigma2 
#>    -16     16 

# Against a numerical derivative of the response Hessian.
f <- function(v) distrib_hess_y(d, y[1], list(mu = v[1], sigma2 = v[2]))
numDeriv::grad(f, c(2, 1))
#> [1] -16  16
```
