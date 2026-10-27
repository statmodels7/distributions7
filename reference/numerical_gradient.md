# Numerical Gradient of the Log-Density

Computes the gradient of the log-density with respect to each parameter
by central finite differences of `distrib_pdf(..., log = TRUE)`. This
powers the default
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
method for distributions that do not implement an analytical gradient:
any `distrib` subclass that defines only `distrib_pdf` receives its
score function from this default at no further cost to its author.

## Usage

``` r
numerical_gradient(distrib, y, theta, h_rel = .Machine$double.eps^(1/3))
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters (each element has length 1 or `length(y)`).

- h_rel:

  Numeric, the relative step size. Defaults to
  `.Machine$double.eps^(1/3)` (optimal for central differences).

## Value

A named list (one element per parameter) of gradient vectors.

## Details

Each component is one central difference of the log-density in its own
parameter,

\$\$l^{(i)} = \frac{\partial}{\partial \theta_i} \log f(y; \theta)
\method \frac{\log f(y; \theta + h_i e_i) - \log f(y; \theta - h_i
e_i)}{2 h_i},\$\$

so the cost is two density evaluations per parameter when the two
candidate steps of
[`fd_steps()`](https://statmodels7.github.io/distributions7/reference/fd_steps.md)
coincide, and more when
[`fd_stable_step()`](https://statmodels7.github.io/distributions7/reference/fd_stable_step.md)
compares them. Truncation is of order \\h^{2}\\ and rounding of order
\\\varepsilon / h\\, which the default \\h \propto \varepsilon^{1/3}\\
balances.

Steps are scaled by `max(1, |theta|)` and kept inside the boundaries of
`distrib@params_bounds`, by
[`fd_stable_step()`](https://statmodels7.github.io/distributions7/reference/fd_stable_step.md),
which reads both of
[`fd_steps()`](https://statmodels7.github.io/distributions7/reference/fd_steps.md)'s
candidates and keeps the one that agrees better with itself at half its
own step. Accuracy is roughly `eps^(2/3)` (about 10 significant digits
at best): sufficient for optimization, but slower and less precise than
an analytical implementation.

## See also

[`numerical_hessian()`](https://statmodels7.github.io/distributions7/reference/numerical_hessian.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)

## Examples

``` r
numerical_gradient(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#> $mu
#> [1] -1  0  1
#> 
#> $sigma
#> [1]  5.500279e-11 -1.000000e+00  5.500279e-11
#> 
```
