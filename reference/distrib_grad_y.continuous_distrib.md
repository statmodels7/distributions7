# Default Response Gradient for Continuous Distributions

Computes \\\partial \ell / \partial y\\ by the two-point central
difference \\\\\ell(y+h) - \ell(y-h)\\/(2h)\\ on
`distrib_pdf(..., log = TRUE)`, through
[`numerical_grad_y()`](https://statmodels7.github.io/distributions7/reference/numerical_grad_y.md).
Registering the fallback on `continuous_distrib` gives a response
gradient to a family that supplies a density and nothing else.

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

The step is \\h = \varepsilon^{1/3}\max(1, \|y\|)\\, shrunk near a
finite bound by
[`fd_steps_y()`](https://statmodels7.github.io/distributions7/reference/fd_steps_y.md)
so that both evaluation points stay inside the support. The stencil's
truncation error is \\O(h^2)\\ and its rounding error
\\O(\varepsilon/h)\\, so that exponent minimizes their sum; measured
against the gamma's closed form it delivers about \\3\times 10^{-10}\\.
The cost is two vectorized
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
calls.

EVERY family the package ships writes its own response gradient, so this
method answers for user-defined families alone. A DISCRETE family has no
derivative in the response, so no default is registered there and the
call raises with the class named.

## Notation

\\\ell\\ is the log-density of one observation, \\y\\ the response,
\\h\\ the finite-difference step and \\\varepsilon\\ the machine
epsilon, `.Machine$double.eps`.

## See also

[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic,
[`numerical_grad_y()`](https://statmodels7.github.io/distributions7/reference/numerical_grad_y.md),
which does the work, and
[`distrib_hess_y.continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.continuous_distrib.md)
for the second order.

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
d <- Toy(distrib_name = "toy", dimension = "univariate",
         bounds = c(-Inf, Inf), params = c("mu", "sigma"),
         params_interpretation = c("location", "scale"), n_params = 2L,
         params_bounds = list(mu = c(-Inf, Inf), sigma = c(0, Inf)),
         link_params = list(mu = linkfunctions7::identity_link(),
                            sigma = linkfunctions7::log_link()))
y <- c(-1, 0.4, 1)
theta <- list(mu = 0.3, sigma = 1.2)

attr(S7::method(distrib_grad_y, Toy), "signature")[[1]]
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
distrib_grad_y(d, y, theta)
#> [1]  1.0138889 -0.8472222 -0.9305556

# Against a Richardson-extrapolated derivative of the same log-density,
# which shares no arithmetic with the stencil.
f <- function(v) distrib_pdf(d, v, theta, log = TRUE)
vapply(y, function(v) numDeriv::grad(f, v), numeric(1))
#> [1]  1.0138889 -0.8472222 -0.9305556

# A discrete family has no method at all.
try(distrib_grad_y(poisson_distrib(), 1:3, list(mu = 2)))
#> Error : Can't find method for `distrib_grad_y(<distributions7::PoissonDistrib>)`.
```
