# Numerical Hyperparameter Hessians of the Response Derivatives

Computes a second-order mixed component by one central difference, in
each parameter, of a first-order component. It supplies the pairs that
have no closed form in the methods for location-scale families with a
shape parameter and for families with a scale and no location. The
default methods for a continuous family use
[`mixed_tensor_derivatives()`](https://statmodels7.github.io/distributions7/reference/mixed_tensor_derivatives.md)
instead.

## Usage

``` r
numerical_theta2_y(distrib, y, theta, inner, h_rel = .Machine$double.eps^(1/3))
```

## Arguments

- distrib:

  A distribution object. Its `params` and `params_bounds` set the
  enumeration and the steps.

- y:

  A numeric vector of observations. The function does not use it,
  because `inner` has already closed over it.

- theta:

  A named list of parameters, the point at which to differentiate.

- inner:

  A function that takes `theta` as its only argument and returns the
  first-order components as a list with one entry per parameter,
  typically `function(th) distrib_cross_y(distrib, y, th)`.

- h_rel:

  The relative step, defaulting to `.Machine$double.eps^(1/3)`, which is
  the optimal exponent for a central first difference.

## Value

A named list with one numeric vector per unordered pair of parameters,
keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## Details

The reference is
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
or
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md),
so a distribution with a closed form for those requires exactly one
difference and never differences a difference.

A mixed pair is differenced in both orders and the two results are
averaged. They agree in exact arithmetic and differ in floating point
because the two steps have different sizes. The average is symmetric in
the pair.

The step is chosen by
[`fd_stable_step()`](https://statmodels7.github.io/distributions7/reference/fd_stable_step.md)
from the two candidates of
[`fd_steps()`](https://statmodels7.github.io/distributions7/reference/fd_steps.md),
so that both evaluation points of the central difference of a parameter
near a bound stay inside its domain.

## Notation

\\\ell\\ is the log-density of one observation, \\y\\ the response,
\\\theta_i\\ a distribution parameter, \\\eta_i\\ its value on the
unconstrained scale and \\h_i = g_i^{-1}\\ the inverse link carrying one
to the other.

## See also

[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
and
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md),
the two generics that it serves, and
[`fd_steps()`](https://statmodels7.github.io/distributions7/reference/fd_steps.md)
for the step rule.

## Examples

``` r
d <- gaussian1_distrib()
y <- c(-1, 0, 2)
theta <- list(mu = 0.4, sigma = 1.3)

num <- numerical_theta2_y(d, y, theta,
                          function(th) distrib_cross_y(d, y, th))
num$sigma_sigma
#> [1]  2.9410735  0.8403067 -3.3612268

# The result is compared with the gaussian's own closed form, for which this
# family does not need the fallback.
max(abs(unlist(num) - unlist(distrib_grad_y_hess(d, y, theta))))
#> [1] 4.191043e-10

# A family that takes the default continuous method is checked against
# numDeriv.
dg <- gamma2_distrib()
yg <- c(0.5, 1, 2)
thg <- list(mu = 2, sigma2 = 1)
g <- distrib_grad_y_hess(dg, yg, thg)
f <- function(v) distrib_grad_y(dg, yg[1], list(mu = v[1], sigma2 = v[2]))
numDeriv::hessian(f, c(2, 1))
#>      [,1] [,2]
#> [1,]    4   -7
#> [2,]   -7   12
c(g$mu_mu[1], g$mu_sigma2[1], g$sigma2_sigma2[1])
#> [1]  4 -7 12
```
