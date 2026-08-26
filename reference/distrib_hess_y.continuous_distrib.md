# Default Response Hessian for Continuous Distributions

Computes \\\partial^2 \ell / \partial y^2\\ by the three-point central
stencil \\\\\ell(y+h) - 2\ell(y) + \ell(y-h)\\/h^2\\ on
`distrib_pdf(..., log = TRUE)`, through
[`numerical_hess_y()`](https://statmodels7.github.io/distributions7/reference/numerical_hess_y.md).
It is one stencil on the log-density and never a difference of
[`distrib_grad_y.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.continuous_distrib.md),
which would be two layers of rounding on the same variable.

## Arguments

- distrib:

  A `continuous_distrib` object with no closed form of its own.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector as long as `y`.

## Details

The step is \\h = \varepsilon^{1/4}\max(1, \|y\|)\\, shrunk near a
finite bound by
[`fd_steps_y()`](https://statmodels7.github.io/distributions7/reference/fd_steps_y.md).
The longer exponent is what the second difference needs: its rounding
error is \\O(\varepsilon/h^2)\\ where the first order's is
\\O(\varepsilon/h)\\, so the optimum sits at a step twenty times longer
and the accuracy reached is coarser, measured at about \\4\times
10^{-7}\\ against the gamma's closed form where the first order reaches
\\3\times 10^{-10}\\. The cost is three vectorized
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
calls.

As at first order, every family the package ships writes its own, and a
DISCRETE family registers no default.

## Notation

\\\ell\\ is the log-density of one observation, \\y\\ the response,
\\h\\ the finite-difference step and \\\varepsilon\\ the machine
epsilon, `.Machine$double.eps`.

## See also

[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic,
[`numerical_hess_y()`](https://statmodels7.github.io/distributions7/reference/numerical_hess_y.md),
which does the work, and
[`distrib_grad_y.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.continuous_distrib.md)
for the first order.

## Examples

``` r
# A family that supplies a density and nothing else, which is what the
# fallback exists for: every family the package ships writes its own
# response derivatives.
Toy <- S7::new_class("Toy", parent = continuous_distrib)
S7::method(distrib_pdf, Toy) <- function(distrib, y, theta, ...) {
  z <- (y - theta$mu) / theta$sigma
  v <- -log(2 * theta$sigma) - abs(z) - 0.1 * z^2
  if (isTRUE(list(...)$log)) v else exp(v)
}
#> Overwriting method distrib_pdf(<Toy>)
d <- Toy(distrib_name = "toy", dimension = "univariate",
         bounds = c(-Inf, Inf), params = c("mu", "sigma"),
         params_interpretation = c("location", "scale"), n_params = 2L,
         params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
         link_params = list(mu = linkfunctions7::identity_link(),
                            sigma = linkfunctions7::log_link()))
y <- c(-1, 0.4, 1)
theta <- list(mu = 0.3, sigma = 1.2)

distrib_hess_y(d, y, theta)
#> [1] -0.1388889 -0.1388889 -0.1388889

# Against a Richardson-extrapolated second derivative of the log-density.
f <- function(v) distrib_pdf(d, v, theta, log = TRUE)
vapply(y, function(v) numDeriv::hessian(f, v)[1, 1], numeric(1))
#> [1] -0.1388889 -0.1388889 -0.1388889

# A discrete family has no method at all.
try(distrib_hess_y(poisson_distrib(), 1:3, list(mu = 2)))
#> Error : Can't find method for `distrib_hess_y(<distributions7::PoissonDistrib>)`.
```
