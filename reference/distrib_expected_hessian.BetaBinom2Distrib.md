# Beta-Binomial Expected Hessian in Its Shapes

Returns the expectation of the observed Hessian under the model,
\\\sum\_{k=0}^{n}
P(Y=k)\\\partial^2\ell/\partial\theta_i\partial\theta_j\\ evaluated at
\\y = k\\. The family being discrete on \\\\0, \dots, n\\\\, that is an
**exact finite sum** of at most \\n+1\\ terms, so the answer is the
expectation to machine precision. `approx` and `nsim` are therefore
ignored: every strategy returns the same three numbers.

The mixed entry does not vanish, so the two shapes are not orthogonal
and their estimates are asymptotically correlated.

## Arguments

- distrib:

  A `BetaBinom2Distrib` object, from
  [`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md).

- y:

  A numeric vector of counts. Only its length is read, the expectation
  not depending on the data; the values are ignored.

- theta:

  A named list with components `alpha` and `beta`, each a numeric vector
  of length 1. Both must be strictly positive. One weighted sum is built
  for the whole call, so a parameter varying by observation is not
  supported here.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored here, the expectation being an exact sum. Accepted so that the
  signature matches the generic's, where it selects between
  `"bartlett"`, `"integrate"`, `"mc"` and `"opg"`.

- nsim:

  Ignored here, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of three numeric vectors, `alpha_alpha`, `beta_beta` and
`alpha_beta`, each of length `length(y)` and each constant along it.

## See also

[`betabinom2_expected()`](https://statmodels7.github.io/distributions7/reference/betabinom2_expected.md)
for the summation,
[`distrib_hessian.BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.BetaBinom2Distrib.md)
for the quantity this is the expectation of, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- betabinom2_distrib(size = 10)
th <- list(alpha = 2, beta = 3)
eh <- distrib_expected_hessian(d, 0:10, th)
vapply(eh, function(v) v[1], numeric(1))
#> alpha_alpha   beta_beta  alpha_beta 
#>  -0.2589858  -0.1111875   0.1523847 

# It is the mass-weighted sum of the observed Hessian over the support,
# written out here by hand and agreeing exactly.
w <- distrib_pdf(d, 0:10, th)
vapply(distrib_hessian(d, 0:10, th), function(v) sum(w * v), numeric(1))
#> alpha_alpha   beta_beta  alpha_beta 
#>  -0.2589858  -0.1111875   0.1523847 

# The strategy argument is inert, the expectation being an exact sum.
identical(eh, distrib_expected_hessian(d, 0:10, th, approx = "mc",
                                       nsim = 50))
#> [1] TRUE
```
