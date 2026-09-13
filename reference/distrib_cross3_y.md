# Mixed Third-Response Parameter Derivatives of the Log-Density

Computes \\\partial^4 \ell / \partial y^3\\ \partial \theta_i\\, one
component per parameter, each a vector along `y`. It is
[`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
differentiated once in each parameter.

## Usage

``` r
distrib_cross3_y(distrib, y, theta, scale = c("parameter", "link"), ...)
```

## Arguments

- distrib:

  A distribution object inheriting from `continuous_distrib`.

- y:

  A numeric vector of observations.

- theta:

  A named list, or named numeric vector, of distribution parameters.
  Each must have length 1 or `length(y)`.

- scale:

  Either `"parameter"` (default) for derivatives with respect to the
  parameters \\\theta\\ on their natural (constrained) scale, or
  `"link"` for derivatives with respect to the unconstrained linear
  predictors \\\eta = g(\theta)\\ defined by `distrib@link_params`. See
  [`link_scale_derivatives()`](https://statmodels7.github.io/distributions7/reference/link_scale_derivatives.md).

- ...:

  Passed to the method.

## Value

A named list with one numeric vector per parameter, each as long as `y`,
keyed by `distrib@params`.

## What consumes it

A penalty \\-\log f(D\beta;\theta)\\ has the Hessian \\S =
-D'\mathrm{diag}(\ell^{(yy)})D\\ in the coefficients, whose derivative
along a direction \\v\\ is \\-D'\mathrm{diag}(\ell^{(yyy)} \odot Dv)D\\.
The second derivative of a marginal criterion in the hyperparameters
reads how that matrix moves with each hyperparameter,
\\-D'\mathrm{diag}(\partial\_{\theta_i}\ell^{(yyy)} \odot Dv)D\\, and
this is the vector placed on its diagonal.

## Where the numbers come from

For a family whose response enters only as \\y - \mu\\ the response
derivative of order three is \\-\partial^3\ell/\partial\mu^3\\, so
\$\$\frac{\partial^4\ell}{\partial y^3\\\partial\theta_i} =
-\frac{\partial^4\ell}{\partial\mu^3\\\partial\theta_i},\$\$ a component
of
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
with a sign. The fourteen location families
[`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
serves by the same identity are served this way. Every other continuous
family takes one central difference of
[`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
in each parameter, through
[`numerical_cross3_y()`](https://statmodels7.github.io/distributions7/reference/numerical_cross3_y.md).

## The link scale

The component for \\\eta_i\\ is the parameter-scale component multiplied
by \\h_i'(\eta_i)\\: the response derivatives are untouched by a
reparametrization of \\\theta\\, so only the first-order diagonal chain
rule enters, as for
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
and
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md).

## Notation

\\\ell\\ is the log-density of one observation, \\y\\ the response,
\\\theta_i\\ a distribution parameter, \\\eta_i\\ its value on the
unconstrained scale and \\h_i = g_i^{-1}\\ the inverse link carrying one
to the other. \\\ell^{(yyy)}\\ is \\\partial^3\ell/\partial y^3\\.

## See also

[`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
for the quantity being differentiated,
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
for the order below, and
[`numerical_cross3_y()`](https://statmodels7.github.io/distributions7/reference/numerical_cross3_y.md)
for the fallback.

## Examples

``` r
d <- student_t1_distrib()
y <- c(-1, 0, 2)
theta <- list(mu = 0.4, sigma = 1.3, nu = 6)
distrib_cross3_y(d, y, theta)
#> $mu
#> [1]  0.02465857 -0.34745651  0.07487536
#> 
#> $sigma
#> [1]  0.7000737  0.4647780 -0.6109113
#> 
#> $nu
#> [1]  0.03085952  0.02819720 -0.02265895
#> 

# Against a numerical derivative of the analytic third response derivative.
f <- function(v) distrib_deriv3_y(d, y[1], list(mu = 0.4, sigma = v, nu = 6))
numDeriv::grad(f, 1.3)
#> [1] 0.7000737

# A family whose response is not a location takes the difference.
distrib_cross3_y(gamma1_distrib(), c(0.5, 1, 2), list(mu = 1.5, phi = 0.8))
#> $mu
#> [1] 0 0 0
#> 
#> $phi
#> [1] -25.000000  -3.125000  -0.390625
#> 
```
