# Mixed Second-Response Parameter Derivatives of the Log-Density

Computes \\\partial^3 \ell / \partial y^2\\ \partial \theta_i\\, one
component per parameter, each a vector along `y`. Where
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
says how curved the log-density is in the RESPONSE, this says how that
curvature moves with each parameter.

## Usage

``` r
distrib_cross2_y(distrib, y, theta, scale = c("parameter", "link"), ...)
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

It completes the surface a marginal likelihood needs. A penalty is
\\-\log f(D\beta;\theta)\\, so its Hessian in the coefficients is
\\-D'\mathrm{diag}(\ell^{(yy)})D\\, and the derivative of that in
\\\theta\\ is this component placed the same way. Differentiating the
Laplace approximation once needs it; differentiating twice needs
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md),
which is this quantity one order further in the parameters.

## The link scale

The component for \\\eta_i\\ is the parameter-scale component multiplied
by \\h_i'(\eta_i)\\. The response derivatives are untouched by a
reparametrization of \\\theta\\, so only the first-order diagonal chain
rule enters, exactly as for
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md).

## Where the numbers come from

A distribution with a closed form provides it. The rest take one central
difference of
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
in each parameter, through
[`numerical_cross2_y()`](https://statmodels7.github.io/distributions7/reference/numerical_cross2_y.md).
The reference is the analytic response Hessian, so a family carrying one
pays for exactly one difference.

Continuous distributions only, as with every response derivative: a
discrete family has no method and the call raises with its class named.

## Notation

\\\ell\\ is the log-density of one observation, \\y\\ the response,
\\\theta_i\\ a distribution parameter, \\\eta_i\\ its value on the
unconstrained scale and \\h_i = g_i^{-1}\\ the inverse link carrying one
to the other. \\\ell^{(yy)}\\ is \\\partial^2\ell/\partial y^2\\.

## See also

[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the quantity being differentiated,
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
for the first-response counterpart,
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md)
for the next order in \\\theta\\, and
[`numerical_cross2_y()`](https://statmodels7.github.io/distributions7/reference/numerical_cross2_y.md)
for the fallback.

## Examples

``` r
d <- gaussian1_distrib()
y <- c(-1, 0, 2)
theta <- list(mu = 0.4, sigma = 1.3)
distrib_cross2_y(d, y, theta)
#> $mu
#> [1] 0 0 0
#> 
#> $sigma
#> [1] 0.9103323 0.9103323 0.9103323
#> 

# The gaussian's response curvature is -1 / sigma^2, which carries no
# location, so the mu component is exactly zero.
c(sigma = 2 / 1.3^3)
#>     sigma 
#> 0.9103323 

# Against a numerical derivative of the response Hessian.
f <- function(v) distrib_hess_y(d, y[1], list(mu = v[1], sigma = v[2]))
numDeriv::grad(f, c(0.4, 1.3))
#> [1] 0.0000000 0.9103323

# On the link scale, sigma riding a log link, the component is multiplied
# by h' = sigma and nothing else.
c(link = distrib_cross2_y(d, y, theta, scale = "link")$sigma[1],
  hand = distrib_cross2_y(d, y, theta)$sigma[1] * 1.3)
#>     link     hand 
#> 1.183432 1.183432 

# A discrete family has no method at all.
try(distrib_cross2_y(poisson_distrib(), 1:3, list(mu = 2)))
#> Error : Can't find method for `distrib_cross2_y(<distributions7::PoissonDistrib>)`.
```
