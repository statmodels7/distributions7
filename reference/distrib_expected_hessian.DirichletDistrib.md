# Dirichlet Expected Hessian

Returns the expectation of the observed Hessian under the model, in
closed form and with no quadrature or simulation:
\$\$\mathbb{E}\[\ell^{(\eta_k\eta_l)}\] = -\phi^2\sum\_{j} t_j
A\_{jk}A\_{jl}, \qquad \mathbb{E}\[\ell^{(\eta_k\phi)}\] =
-\phi\sum\_{j} t_j \mu_j A\_{jk}, \qquad \mathbb{E}\[\ell^{(\phi\phi)}\]
= \psi'(\phi) - \sum\_{j} t_j \mu_j^2,\$\$ with \\t_j =
\psi'(\alpha_j)\\ and \\A = \partial\mu/\partial\eta\\. `approx` and
`nsim` are ignored: every strategy returns the same matrix.

## Arguments

- distrib:

  A `DirichletDistrib` object, from
  [`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md).

- y:

  A numeric matrix of observations. Only its row count is read through
  [`n_obs()`](https://statmodels7.github.io/distributions7/reference/n_obs.md),
  the expectation not depending on the data; the values are ignored.

- theta:

  A named list of parameters on the parameter scale: the mean's free
  values followed by `phi`, each of length 1.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored here, the expectation being exact. Accepted so that the
  signature matches the generic's, where it selects between
  `"bartlett"`, `"integrate"`, `"mc"` and `"opg"`.

- nsim:

  Ignored here, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors, one per distinct entry of the symmetric
matrix and named as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
names them, each of length `n_obs(distrib, y)` and constant along it.

## Why the data-carrying terms drop out

Under the model \\\mathbb{E}\[\log y_j\] = \psi(\alpha_j) -
\psi(\phi)\\, so \\\mathbb{E}\[g_j\] = -\psi(\phi)\\, **the same
constant for every** \\j\\. Each term of the observed Hessian that
carries the data is that constant times a sum over \\j\\ of a column of
\\A\\ or of a second-derivative vector of the simplex, and both of those
sum to zero: differentiating \\\sum_j \mu_j = 1\\ once gives the first,
twice the second. Every such term is therefore exactly zero under
expectation, which is why a family with no location and scale to
separate still has a closed-form information matrix.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean vector,
\\\phi \> 0\\ the concentration, \\\alpha = \phi\mu\\ the shapes,
\\\eta\\ the free vector of the simplex chart, and \\\psi\\, \\\psi'\\
the digamma and trigamma functions.

## See also

[`distrib_hessian.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.DirichletDistrib.md)
for the quantity this is the expectation of,
[`distrib_gradient.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.DirichletDistrib.md)
for the score, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- dirichlet_distrib(3)
th <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 12)
set.seed(1)
Y <- distrib_rng(d, 4, th)
eh <- distrib_expected_hessian(d, Y, th)
vapply(eh, function(v) v[1], numeric(1))
#> mean_alr1_mean_alr1 mean_alr2_mean_alr2             phi_phi mean_alr1_mean_alr2 
#>        -3.308634934        -2.675298951        -0.007740407         1.514683277 
#>       mean_alr1_phi       mean_alr2_phi 
#>         0.013449776        -0.010993516 

# Averaging the observed Hessian over draws reaches the same matrix.
set.seed(3)
Z <- distrib_rng(d, 3e5, th)
vapply(distrib_hessian(d, Z, th), mean, numeric(1))
#> mean_alr1_mean_alr1 mean_alr2_mean_alr2             phi_phi mean_alr1_mean_alr2 
#>        -3.308487877        -2.675658807        -0.007740407         1.514743742 
#>       mean_alr1_phi       mean_alr2_phi 
#>         0.013532593        -0.011055574 

# The two zero sums the derivation rests on, read off the simplex itself.
eta <- c(0.3, -0.2)
colSums(do.call(cbind, parameters7::param_d1(d@param, eta)))
#>         alr1         alr2 
#> 2.775558e-17 2.775558e-17 
vapply(parameters7::param_d2(d@param, eta), sum, numeric(1))
#>     alr1:alr1     alr2:alr2     alr1:alr2 
#>  6.938894e-18  1.387779e-17 -2.081668e-17 

# The strategy argument is inert, the expectation being exact.
identical(eh, distrib_expected_hessian(d, Y, th, approx = "mc", nsim = 50))
#> [1] TRUE
```
