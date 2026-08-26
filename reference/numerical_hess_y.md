# Numerical Second Response Derivative of the Log-Density

Computes \\\partial^2 \ell / \partial y^2\\ by the three-point central
stencil \$\$\frac{\partial^2\ell}{\partial y^2} \approx
\frac{\ell(y+h) - 2\ell(y) + \ell(y-h)}{h^2},\$\$ evaluated on
`distrib_pdf(..., log = TRUE)`. It is what the default
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
method runs for a continuous family with no closed form.

## Usage

``` r
numerical_hess_y(distrib, y, theta, h_rel = .Machine$double.eps^(1/4))
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

  The relative step, defaulting to `.Machine$double.eps^(1/4)`, the
  exponent that minimizes the total error of a central second
  difference.

## Value

A numeric vector as long as `y`.

## Why the step is longer and the answer coarser

The truncation error is again \\O(h^2)\\, but the rounding error is
\\O(\varepsilon/h^2)\\ where the first order's is \\O(\varepsilon/h)\\,
the numerator being a difference of quantities of the same size divided
by \\h^2\\. The total is smallest at \\h \sim \varepsilon^{1/4}\\, which
is the default and is twenty times the first order's step, and the
accuracy attainable is correspondingly coarser. Measured against the
gamma's own closed form at `h_rel` from \\10^{-1}\\ to \\10^{-8}\\:
2.5e-01, 2.4e-03, 2.4e-05, 2.8e-07, 5.4e-06, 6.2e-04, 3.5e-02 and
3.1e+00. The rise past the optimum is far steeper than the first
order's, which is the \\h^{-2}\\ at work.

## The cost

Three evaluations of the log-density, one of them at `y` itself, all
vectorized over `y`.

## Notation

\\\ell\\ is the log-density of one observation, \\y\\ the response,
\\h\\ the finite-difference step and \\\varepsilon\\ the machine
epsilon, `.Machine$double.eps`.

## See also

[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic,
[`numerical_grad_y()`](https://statmodels7.github.io/distributions7/reference/numerical_grad_y.md)
for the first order, and
[`fd_steps_y()`](https://statmodels7.github.io/distributions7/reference/fd_steps_y.md)
for the step rule.

## Examples

``` r
d <- gamma2_distrib()
y <- c(0.5, 1, 2)
theta <- list(mu = 2, sigma2 = 1)

numerical_hess_y(d, y, theta)
#> [1] -12.00  -3.00  -0.75
distrib_hess_y(d, y, theta)
#> [1] -12.00  -3.00  -0.75

# The second difference is the cruder of the two, by about three orders of
# magnitude on the same family and the same points.
c(order2 = max(abs(numerical_hess_y(d, y, theta) -
                   distrib_hess_y(d, y, theta))),
  order1 = max(abs(numerical_grad_y(d, y, theta) -
                   distrib_grad_y(d, y, theta))))
#>       order2       order1 
#> 3.576279e-07 3.142242e-10 

# And it degrades far faster below its optimal step.
err <- function(p) max(abs(numerical_hess_y(d, y, theta, h_rel = 10^(-p)) -
                           distrib_hess_y(d, y, theta)))
setNames(vapply(3:8, err, numeric(1)), paste0("1e-", 3:8))
#>         1e-3         1e-4         1e-5         1e-6         1e-7         1e-8 
#> 2.400026e-05 2.823417e-07 5.433777e-06 6.227178e-04 3.481759e-02 3.118216e+00 
```
