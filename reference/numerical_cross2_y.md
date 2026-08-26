# Numerical Mixed Second-Response Parameter Derivatives

Computes \\\partial^3 \ell / \partial y^2\\ \partial \theta_i\\ by one
central difference of
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
in each parameter. It is what the default
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
method runs for a continuous family with no closed form of its own.

## Usage

``` r
numerical_cross2_y(
  distrib,
  y,
  theta,
  h_rel = .Machine$double.eps^(1/3),
  which = NULL
)
```

## Arguments

- distrib:

  A distribution object inheriting from `continuous_distrib`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, the point to differentiate at.

- h_rel:

  The relative step, defaulting to `.Machine$double.eps^(1/3)`, the
  optimal exponent for a central first difference.

- which:

  An optional character vector naming a subset of `distrib@params`. The
  default, `NULL`, computes every component; a subset returns only
  those, keyed the same way, and costs one difference each.

## Value

A named list with one numeric vector per requested parameter, each as
long as `y`.

## Details

The reference is the response HESSIAN, one order up from the
log-density, so a distribution carrying an analytic `distrib_hess_y`
pays for exactly one finite-difference layer. Where that Hessian is
itself a fallback the two differences act on different variables and
compose into a single mixed stencil, rather than into the nested
differencing of one variable the package forbids.

The step is
[`fd_steps()`](https://statmodels7.github.io/distributions7/reference/fd_steps.md)'s,
so a parameter near a bound is differenced inward instead of across it.

## Notation

\\\ell\\ is the log-density of one observation, \\y\\ the response,
\\\theta_i\\ a distribution parameter, \\\eta_i\\ its value on the
unconstrained scale and \\h_i = g_i^{-1}\\ the inverse link carrying one
to the other. \\\ell^{(yy)}\\ is \\\partial^2\ell/\partial y^2\\.

## See also

[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md),
the generic it serves,
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the quantity it differences, and
[`fd_steps()`](https://statmodels7.github.io/distributions7/reference/fd_steps.md)
for the step rule.

## Examples

``` r
d <- gaussian1_distrib()
y <- c(-1, 0, 2)
theta <- list(mu = 0.4, sigma = 1.3)

num <- numerical_cross2_y(d, y, theta)
num$sigma
#> [1] 0.9103323 0.9103323 0.9103323

# Against the gaussian's own closed form, which it does not need this for.
max(abs(unlist(num) - unlist(distrib_cross2_y(d, y, theta))))
#> [1] 6.796563e-11

# One parameter at a time, when only one is wanted.
numerical_cross2_y(d, y, theta, which = "sigma")
#> $sigma
#> [1] 0.9103323 0.9103323 0.9103323
#> 

# A family that does take this route, checked against numDeriv.
dg <- gamma2_distrib()
yg <- c(0.5, 1, 2)
vapply(distrib_cross2_y(dg, yg, list(mu = 2, sigma2 = 1)),
       function(z) z[1], numeric(1))
#>     mu sigma2 
#>    -16     16 
f <- function(v) distrib_hess_y(dg, yg[1], list(mu = v[1], sigma2 = v[2]))
numDeriv::grad(f, c(2, 1))
#> [1] -16  16
```
