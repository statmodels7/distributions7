# Dirichlet Observed Hessian

Computes the distinct second derivatives of the Dirichlet log-density
with respect to the mean's free values and the concentration, one value
per observation, in closed form. With \\g_j = \log y_j -
\psi(\alpha_j)\\, \\t_j = \psi'(\alpha_j)\\, \\A =
\partial\mu/\partial\eta\\ and \\B\_{\cdot,kl} =
\partial^2\mu/\partial\eta_k\partial\eta_l\\, \$\$\ell^{(\eta_k\eta_l)}
= \phi\sum\_{j}\left(-t_j\phi A\_{jk}A\_{jl} + g_j B\_{j,kl}\right),
\qquad \ell^{(\eta_k\phi)} = -\phi\sum\_{j} t_j \mu_j A\_{jk} +
\sum\_{j} g_j A\_{jk},\$\$ \$\$\ell^{(\phi\phi)} = \psi'(\phi) -
\sum\_{j} t_j \mu_j^2.\$\$ The pure-concentration entry is **free of the
data**, the family being an exponential family in \\\log y\\ whose
concentration multiplies a statistic linearly; it therefore equals its
own expectation at every observation. The other two carry the data
through \\g\\.

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

A named list of numeric vectors, one per distinct entry of the symmetric
matrix and named as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
names them, each of length `nrow(y)`. For \\p\\ coordinates there are
\\p(p+1)/2\\ of them.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean vector,
\\\phi \> 0\\ the concentration, \\\alpha = \phi\mu\\ the shapes,
\\\eta\\ the free vector of the simplex chart, and \\\psi\\, \\\psi'\\
the digamma and trigamma functions.

## See also

[`distrib_gradient.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.DirichletDistrib.md)
for the score,
[`distrib_expected_hessian.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.DirichletDistrib.md)
for the expectation of this quantity,
[`distrib_deriv3.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.DirichletDistrib.md)
for the order above, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- dirichlet_distrib(3)
th <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 12)
set.seed(1)
Y <- distrib_rng(d, 4, th)
h <- distrib_hessian(d, Y, th)
names(h)
#> [1] "mean_alr1_mean_alr1" "mean_alr2_mean_alr2" "phi_phi"            
#> [4] "mean_alr1_mean_alr2" "mean_alr1_phi"       "mean_alr2_phi"      

# The pure-concentration entry is the same number at every observation.
h$phi_phi
#> [1] -0.007740407 -0.007740407 -0.007740407 -0.007740407

# numDeriv on the summed log-density reproduces the summed matrix.
fn <- function(v)
  sum(distrib_pdf(d, Y, as.list(stats::setNames(v, d@params)), log = TRUE))
H <- numDeriv::hessian(fn, unlist(th))
rbind(numeric = c(H[1, 1], H[2, 2], H[3, 3], H[1, 2], H[1, 3], H[2, 3]),
      closed = vapply(h, sum, numeric(1)))
#>         mean_alr1_mean_alr1 mean_alr2_mean_alr2     phi_phi mean_alr1_mean_alr2
#> numeric           -12.64827           -10.24125 -0.03096163            4.629503
#> closed            -12.64827           -10.24125 -0.03096163            4.629503
#>         mean_alr1_phi mean_alr2_phi
#> numeric     0.3839644    0.03534577
#> closed      0.3839644    0.03534577
```
