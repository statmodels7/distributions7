# Numerical Gradient of the Log-Density with Respect to the Response

Computes \\\partial \ell / \partial y\\ by the two-point central
difference \$\$\frac{\partial\ell}{\partial y} \approx \frac{\ell(y+h) -
\ell(y-h)}{2h},\$\$ evaluated on `distrib_pdf(..., log = TRUE)`. It is
what the default
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
method runs for a continuous family with no closed form.

## Usage

``` r
numerical_grad_y(distrib, y, theta, h_rel = .Machine$double.eps^(1/3))
```

## Arguments

- distrib:

  An object inheriting from `continuous_distrib`. A discrete family has
  no derivative in the response and registers no method.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- h_rel:

  The relative step, defaulting to `.Machine$double.eps^(1/3)`, the
  exponent that minimizes the total error of a central first difference.

## Value

A numeric vector as long as `y`.

## The step and the accuracy it buys

The stencil's truncation error is \\O(h^2)\\ and its rounding error
\\O(\varepsilon/h)\\, so the total is smallest at \\h \sim
\varepsilon^{1/3}\\, which is the default. Measured against the gamma's
own closed form at `h_rel` from \\10^{-1}\\ to \\10^{-8}\\, the error is
8.2e-02, 8.0e-04, 8.0e-06, 8.0e-08, 8.0e-10, 1.4e-10, 8.2e-10 and
2.0e-08: it divides by a hundred per decade while the truncation
dominates, turns near the default and rises again as the rounding takes
over. A step ten times shorter than the default is therefore worse, not
better.

The step itself comes from
[`fd_steps_y()`](https://statmodels7.github.io/distributions7/reference/fd_steps_y.md),
which scales it by \\\max(1, \|y\|)\\ and shrinks it near a finite
bound.

## The cost

Two evaluations of the log-density, both vectorized over `y`, so one
call costs two
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
calls whatever the sample size. A family with a closed form pays
neither.

## Notation

\\\ell\\ is the log-density of one observation, \\y\\ the response,
\\h\\ the finite-difference step and \\\varepsilon\\ the machine
epsilon, `.Machine$double.eps`.

## See also

[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
for the generic,
[`numerical_hess_y()`](https://statmodels7.github.io/distributions7/reference/numerical_hess_y.md)
for the second order, and
[`fd_steps_y()`](https://statmodels7.github.io/distributions7/reference/fd_steps_y.md)
for the step rule.

## Examples

``` r
d <- gamma2_distrib()
y <- c(0.5, 1, 2)
theta <- list(mu = 2, sigma2 = 1)

numerical_grad_y(d, y, theta)
#> [1]  4.0  1.0 -0.5

# Against the gamma's own closed form: the difference is good to about
# 3e-10, not to machine precision.
distrib_grad_y(d, y, theta)
#> [1]  4.0  1.0 -0.5
max(abs(numerical_grad_y(d, y, theta) - distrib_grad_y(d, y, theta)))
#> [1] 3.142242e-10

# Too short a step is worse than too long. The error falls as h^2 and then
# rises as eps / h, and the default sits at the turn.
err <- function(p) max(abs(numerical_grad_y(d, y, theta, h_rel = 10^(-p)) -
                           distrib_grad_y(d, y, theta)))
setNames(vapply(1:8, err, numeric(1)), paste0("1e-", 1:8))
#>         1e-1         1e-2         1e-3         1e-4         1e-5         1e-6 
#> 8.197662e-02 8.001921e-04 8.000019e-06 7.999779e-08 8.033609e-10 1.397780e-10 
#>         1e-7         1e-8 
#> 8.182894e-10 2.009904e-08 
```
