# A Second-Order Mixed Derivative Through a Map

Carries the parent's paired components onto a new parametrization by the
second-order chain rule. It is the arithmetic behind every
[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
and
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md)
method of a reparametrized family, and it differentiates nothing itself.

## Usage

``` r
mapped_theta2(distrib, parent, th_par, maps, y, first, second)
```

## Arguments

- distrib:

  The distribution in the NEW parametrization. Only its `params` are
  read, to key the result.

- parent:

  The parent distribution, whose `params` key the inputs.

- th_par:

  A named list of the parent's parameters, evaluated at the new ones.

- maps:

  The map's keyed partial tables, from
  [`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md):
  one entry per parent parameter, each a list keyed by new-parameter
  index.

- y:

  A numeric vector of observations. Used only for its length, through
  the components.

- first:

  The parent's first-order components, keyed by the parent's parameters.

- second:

  The parent's second-order components, keyed by the parent's parameter
  pairs. Either ordering of a pair's key is accepted.

## Value

A named list with one numeric vector per unordered pair of the NEW
parameters, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## The expansion

A derivative in the RESPONSE does not see the parameters, so a
reparametrization acts on the \\\theta\\ side alone and the ordinary
two-term expansion applies: \$\$\frac{\partial^2 f}{\partial\alpha_a
\partial\alpha_b} = \sum\_{i,j}
\frac{\partial\theta_i}{\partial\alpha_a}
\frac{\partial\theta_j}{\partial\alpha_b} f\_{ij} + \sum_i
\frac{\partial^2\theta_i}{\partial\alpha_a\partial\alpha_b} f_i,\$\$
with \\f_i\\ the parent's first-order mixed component and \\f\_{ij}\\
its second-order one. Both partial tables are the MAP's, already keyed
by
[`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md)
under the indices of the new parameters, `"1"` for a first partial and
`"1,2"` for a second, so nothing is differentiated here.

## What `first` and `second` are

They are the same quantity at two orders in \\\theta\\, and which
quantity depends on the caller:
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
with
[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
for the third-order derivative, and
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
with
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md)
for the fourth.

A missing key in a map's table is an exact zero and is skipped, so a map
that is affine in one of its coordinates costs nothing for it.

## Notation

\\\ell\\ is the log-density of one observation, \\\theta_i\\ a parameter
of the PARENT and \\\alpha_a\\ one of the new parametrization, so that
\\\theta = \theta(\alpha)\\ is the map. \\\ell^{(y)}\\ and
\\\ell^{(yy)}\\ are the first and second derivatives in the response.

## See also

[`mapped_cross2_y()`](https://statmodels7.github.io/distributions7/reference/mapped_cross2_y.md)
for the first-order counterpart,
[`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md)
for the partial tables, and
[`distrib_grad_y_hess.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.ReparamContinuousDistrib.md),
its caller.

## Examples

``` r
# gaussian2 carries (mu, sigma2) where its parent carries (mu, sigma).
d <- gaussian2_distrib()
y <- c(-0.7, 0.3, 1.4)
theta <- list(mu = 0.3, sigma2 = 1.44)

g <- distrib_grad_y_hess(d, y, theta)
vapply(g, function(z) z[1], numeric(1))
#>         mu_mu sigma2_sigma2     mu_sigma2 
#>     0.0000000     0.6697960    -0.4822531 

# Against a numerical Hessian of the response gradient in the NEW
# parameters, which shares none of the map's arithmetic.
f <- function(v) distrib_grad_y(d, y[1], list(mu = v[1], sigma2 = v[2]))
numDeriv::hessian(f, c(0.3, 1.44))
#>               [,1]       [,2]
#> [1,]  9.210739e-13 -0.4822531
#> [2,] -4.822531e-01  0.6697960
```
