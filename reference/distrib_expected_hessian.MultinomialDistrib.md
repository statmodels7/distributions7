# Multinomial Expected Hessian

Returns the expectation of the observed Hessian under the model, in
closed form and with no quadrature, simulation or summation over the
support: \$\$\mathbb{E}\[\ell^{(\eta_k\eta_l)}\] = -n\sum\_{j=1}^{p}
\dfrac{A\_{jk}A\_{jl}}{p_j},\$\$ with \\A = \partial p/\partial\eta\\.
`approx` and `nsim` are ignored: every strategy returns the same matrix.

## Arguments

- distrib:

  A `MultinomialDistrib` object, from
  [`multinomial_distrib()`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md).

- y:

  A numeric matrix of observations. Only its row count is read through
  [`n_obs()`](https://statmodels7.github.io/distributions7/reference/n_obs.md),
  the expectation not depending on the data; the values are ignored.

- theta:

  A named list of the simplex's free values on the parameter scale, each
  of length 1.

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

## Why the second-derivative term vanishes

Under the model \\\mathbb{E}\[y_j\] = n p_j\\, so the term of the
observed Hessian carrying \\B\\ becomes \\n\sum_j B\_{j,kl}\\. The
probabilities sum to one at every free vector, so every derivative of
that sum is zero and the whole term drops out, leaving the display
above. The same argument, one order lower, gives the score mean zero.

## The information is the covariance of the first p-1 counts

On the default chart,
[`parameters7::simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.html)'s
additive log-ratio coordinates, the Jacobian is \\A\_{jk} =
p_j(\delta\_{jk} - p_k)\\, so \$\$\sum_j \dfrac{A\_{jk}A\_{jl}}{p_j} =
\sum_j p_j(\delta\_{jk}-p_k)(\delta\_{jl}-p_l) = \delta\_{kl}p_k - p_k
p_l,\$\$ and the expected information
\\-\mathbb{E}\[\ell^{(\eta_k\eta_l)}\]\\ is exactly
\\\operatorname{Cov}(Y_k, Y_l)\\, the leading \\(p-1)\times(p-1)\\ block
of
[`mv_sigma.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MultinomialDistrib.md).
Measured, the two agree to `2e-16` at every dimension and trial count
tried. The identity belongs to that chart; a different simplex
parametrization has a different Jacobian and no such coincidence.

## Notation

\\\ell\\ is the log-mass of one observation, \\p\\ the probability
vector, \\n\\ the trial count, \\\eta\\ the free vector of the simplex
chart and \\A = \partial p/\partial\eta\\ its Jacobian.

## See also

[`distrib_hessian.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.MultinomialDistrib.md)
for the quantity this is the expectation of,
[`mv_sigma.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MultinomialDistrib.md)
for the covariance it coincides with,
[`mv_support.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_support.MultinomialDistrib.md)
for the exact sum it can be checked against, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- multinomial_distrib(3, size = 5)
th <- list(probs_alr1 = 0.3, probs_alr2 = -0.2)
eh <- distrib_expected_hessian(d, matrix(0, 1, 3), th)
vapply(eh, function(v) v[1], numeric(1))
#> probs_alr1_probs_alr1 probs_alr2_probs_alr2 probs_alr1_probs_alr2 
#>            -1.2226293            -0.9581222             0.5503861 

# The exact sum over the support gives the same three numbers.
supp <- mv_support(d, th)
w <- distrib_pdf(d, supp, th)
vapply(distrib_hessian(d, supp, th), function(v) sum(w * v), numeric(1))
#> probs_alr1_probs_alr1 probs_alr2_probs_alr2 probs_alr1_probs_alr2 
#>            -1.2226293            -0.9581222             0.5503861 

# On the default chart the information is the covariance of the first
# p - 1 counts.
S <- mv_sigma(d, th)
rbind(information = -vapply(eh, function(v) v[1], numeric(1)),
      covariance = c(S[1, 1], S[2, 2], S[1, 2]))
#>             probs_alr1_probs_alr1 probs_alr2_probs_alr2 probs_alr1_probs_alr2
#> information              1.222629             0.9581222            -0.5503861
#> covariance               1.222629             0.9581222            -0.5503861

# The strategy argument is inert, the expectation being exact.
identical(eh, distrib_expected_hessian(d, matrix(0, 1, 3), th,
                                       approx = "mc", nsim = 50))
#> [1] TRUE
```
