# The Link Scale of a Second-Order Parameter Derivative

Carries a quantity keyed by parameter pair from the parameter scale onto
the unconstrained one. Both
[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
and
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md)
call it after dispatch, so a method returns the parameter scale and
never has to know about links.

## Usage

``` r
theta2_link_scale(distrib, y, theta, second, first)
```

## Arguments

- distrib:

  A distribution object. Its `link_params` supply \\h'\\ and \\h''\\
  through
  [`inverse_link_derivs()`](https://statmodels7.github.io/distributions7/reference/inverse_link_derivs.md).

- y:

  The response. Unused, and present so that the signature reads like its
  callers'.

- theta:

  A named list of parameters, on the parameter scale.

- second:

  The parameter-scale second-order components, keyed by parameter pair.

- first:

  The parameter-scale first-order components, keyed by parameter:
  [`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
  for the third-order quantity and
  [`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
  for the fourth.

## Value

A named list keyed exactly as `second`.

## Details

The chain rule is DIAGONAL, each parameter carrying its own link, so a
pair \\(i, j)\\ becomes \$\$\frac{\partial^2}{\partial\eta_i
\partial\eta_j} = h_i'(\eta_i)\\ h_j'(\eta_j)
\frac{\partial^2}{\partial\theta_i \partial\theta_j} + \delta\_{ij}\\
h_i''(\eta_i) \frac{\partial}{\partial\theta_i}.\$\$ The response
derivatives do not enter: a reparametrization of \\\theta\\ leaves them
alone, which is why `first` and `second` are the SAME quantity at two
orders in \\\theta\\ and no third argument is needed.

## Notation

\\\ell\\ is the log-density of one observation, \\y\\ the response,
\\\theta_i\\ a distribution parameter, \\\eta_i\\ its value on the
unconstrained scale and \\h_i = g_i^{-1}\\ the inverse link carrying one
to the other.

## See also

[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
and
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md),
its two callers, and
[`inverse_link_derivs()`](https://statmodels7.github.io/distributions7/reference/inverse_link_derivs.md)
for the link derivatives.

## Examples

``` r
d <- gaussian1_distrib()
y <- c(-1, 0, 2)
theta <- list(mu = 0.4, sigma = 1.3)

# What the generic does when scale = "link" is asked for.
second <- distrib_grad_y_hess(d, y, theta)
first <- distrib_cross_y(d, y, theta)
linked <- distributions7:::theta2_link_scale(d, y, theta, second, first)
identical(linked, distrib_grad_y_hess(d, y, theta, scale = "link"))
#> [1] TRUE

# sigma rides a log link, so h' = h'' = sigma and the diagonal pair gains
# the first-order term while the off-diagonal one does not.
c(diagonal = second$sigma_sigma[1] * 1.3^2 + first$sigma[1] * 1.3,
  off = second$mu_sigma[1] * 1 * 1.3)
#>  diagonal       off 
#>  3.313609 -1.183432 
c(linked$sigma_sigma[1], linked$mu_sigma[1])
#> [1]  3.313609 -1.183432
```
