# Dirichlet Score

Computes the first derivatives of the Dirichlet log-density with respect
to the mean's free values and the concentration, one value per
observation, in closed form. Writing \\g_j = \log y_j - \psi(\alpha_j)\\
and \\A = \partial\mu/\partial\eta\\ for the simplex's own Jacobian,
\$\$\dfrac{\partial\ell}{\partial\eta_k} = \phi \sum\_{j} g_j A\_{jk},
\qquad \dfrac{\partial\ell}{\partial\phi} = \psi(\phi) + \sum\_{j} g_j
\mu_j.\$\$ The data enter only through \\\log y\\, the family being an
exponential family in that statistic, so \\g\\ carries the whole of it.

With `scale = "link"` the generic applies the chain rule for the links
the family carries before returning. Every mean coordinate rides the
identity, its free value being unconstrained already, so only the
concentration's component changes.

## Arguments

- distrib:

  A `DirichletDistrib` object, from
  [`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md).

- y:

  A numeric matrix with one row per observation and \\p\\ columns, each
  row on the open simplex. A single observation may be given as a plain
  vector.

- theta:

  A named list of parameters on the parameter scale: the mean's free
  values followed by `phi`, each of length 1.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with one numeric vector per parameter, in the order
`distrib@params` gives, each of length `nrow(y)`.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean vector,
\\\phi \> 0\\ the concentration, \\\alpha = \phi\mu\\ the shapes,
\\\eta\\ the free vector of the simplex chart and \\\psi\\ the digamma
function.

## See also

[`distrib_hessian.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.DirichletDistrib.md)
for the second derivatives,
[`distrib_expected_hessian.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.DirichletDistrib.md)
for their expectation,
[`dir_parts()`](https://statmodels7.github.io/distributions7/reference/dir_parts.md)
for the shared pieces, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- dirichlet_distrib(3)
th <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 12)
set.seed(1)
Y <- distrib_rng(d, 4, th)
g <- distrib_gradient(d, Y, th)
names(g)
#> [1] "mean_alr1" "mean_alr2" "phi"      

# It is the derivative of the log-density, so numDeriv on the summed
# log-density reproduces the summed score.
fn <- function(v)
  sum(distrib_pdf(d, Y, as.list(stats::setNames(v, d@params)), log = TRUE))
rbind(numeric = numDeriv::grad(fn, unlist(th)),
      closed = vapply(g, sum, numeric(1)))
#>         mean_alr1 mean_alr2         phi
#> numeric  3.961983  0.951838 0.004907997
#> closed   3.961983  0.951838 0.004907997

# The summed score vanishes at the maximum likelihood estimate.
set.seed(9)
Z <- distrib_rng(d, 800, th)
mle <- as.list(coef(fit_distrib(d, Z)))
vapply(distrib_gradient(d, Z, mle), sum, numeric(1))
#>     mean_alr1     mean_alr2           phi 
#> -4.977125e-05  4.248382e-05  2.472542e-05 
```
