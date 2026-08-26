# Numerical Hyperparameter Hessians of the Response Derivatives

Computes a second-order mixed component by one central difference, in
each parameter, of an ANALYTIC first-order component. It is the fallback
behind
[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
and
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md)
for a family that provides no closed form.

## Usage

``` r
numerical_theta2_y(distrib, y, theta, inner, h_rel = .Machine$double.eps^(1/3))
```

## Arguments

- distrib:

  A distribution object. Its `params` and `params_bounds` set the
  enumeration and the steps.

- y:

  A numeric vector of observations. Passed to nothing here; `inner` has
  already closed over it.

- theta:

  A named list of parameters, the point to differentiate at.

- inner:

  A function of `theta` alone returning the first-order components as a
  list with one entry per parameter, typically
  `function(th) distrib_cross_y(distrib, y, th)`.

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
so a distribution with a closed form for those pays for exactly one
difference and never differences a difference.

A mixed pair is differenced BOTH WAYS and averaged. The two agree in
exact arithmetic, and in floating point they differ, the two steps being
different sizes, while a second derivative of a scalar has to come out
symmetric. On a gamma the two orders differ by about \\3\times
10^{-10}\\, which is the size of the answer's own error.

The step is
[`fd_steps()`](https://statmodels7.github.io/distributions7/reference/fd_steps.md)'s,
so a parameter near a bound is differenced inward rather than across it.

## Notation

\\\ell\\ is the log-density of one observation, \\y\\ the response,
\\\theta_i\\ a distribution parameter, \\\eta_i\\ its value on the
unconstrained scale and \\h_i = g_i^{-1}\\ the inverse link carrying one
to the other.

## See also

[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
and
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md),
the two generics it serves, and
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

# Against the gaussian's own closed form, which this family does not need
# the fallback for.
max(abs(unlist(num) - unlist(distrib_grad_y_hess(d, y, theta))))
#> [1] 4.191043e-10

# A family that does take the fallback, checked against numDeriv.
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
