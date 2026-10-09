# Default Hyperparameter Hessian of the Response Curvature

Falls back to one tensor stencil on the highest quantity the family
implements itself, through
[`mixed_tensor_derivatives()`](https://statmodels7.github.io/distributions7/reference/mixed_tensor_derivatives.md).
It is
[`distrib_grad_y_hess.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.continuous_distrib.md)
read one order higher in the response, and it makes the fourth-order
mixed derivative available for every continuous family. For
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md),
which has an analytic Hessian in the parameters and no response
derivative, it is a second difference in the response of that Hessian;
before 0.76.0 it differenced a difference in the same parameter, and on
a log-gamma it was wrong by a factor of seven.

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

## See also

[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md)
for the generic,
[`distrib_grad_y_hess.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.continuous_distrib.md)
for the third-order twin, and
[`numerical_theta2_y()`](https://statmodels7.github.io/distributions7/reference/numerical_theta2_y.md),
which does the work.

## Examples

``` r
d <- gamma2_distrib()
y <- c(0.5, 1, 2)
theta <- list(mu = 2, sigma2 = 1)
h <- distrib_hess_y_hess(d, y, theta)
c(h$mu_mu[1], h$mu_sigma2[1], h$sigma2_sigma2[1])
#> [1]  -8  16 -32

# Against a numerical Hessian of the response curvature.
f <- function(v) distrib_hess_y(d, y[1], list(mu = v[1], sigma2 = v[2]))
numDeriv::hessian(f, c(2, 1))
#>      [,1] [,2]
#> [1,]   -8   16
#> [2,]   16  -32
```
