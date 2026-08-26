# Beta-Binomial Expected Hessian

Returns the expectation of the observed Hessian under the model,
\\\sum\_{k=0}^{n} P(Y=k)\\
\partial^2\ell/\partial\theta_i\partial\theta_j\\ evaluated at \\y =
k\\. The support being finite, that is an **exact finite sum** of
\\n+1\\ terms rather than a quadrature or a sample. The answer is the
expectation to machine precision, and a bounded support is what makes
that available.

`approx` and `nsim` are therefore ignored; every strategy returns the
same three numbers. The arithmetic runs in a compiled kernel decomposed
over the elements of the output, so the result does not depend on the
thread count.

## Arguments

- distrib:

  A `BetaBinom1Distrib` object, from
  [`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md).

- y:

  A numeric vector of counts. Only its length is read, the expectation
  not depending on the data; the values are ignored.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `mu` must lie in \\(0, 1)\\ and `sigma` be strictly
  positive.

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

- threads:

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count says. Defaults to `1L`.

## Value

A named list of three numeric vectors, `mu_mu`, `mu_sigma` and
`sigma_sigma`, each of length `length(y)` and each constant along it.

## See also

[`distrib_hessian.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.BetaBinom1Distrib.md)
for the quantity this is the expectation of,
[`distrib_gradient.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.BetaBinom1Distrib.md)
for the score, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- betabinom1_distrib(size = 10)
th <- list(mu = 0.3, sigma = 0.5)
eh <- distrib_expected_hessian(d, 0:10, th)
vapply(eh, function(v) v[1], numeric(1))
#>       mu_mu    mu_sigma sigma_sigma 
#> -12.9931856   0.9716737  -1.2084218 

# It is the mass-weighted sum of the observed Hessian over the support,
# written out here by hand and agreeing exactly.
w <- distrib_pdf(d, 0:10, th)
vapply(distrib_hessian(d, 0:10, th), function(v) sum(w * v), numeric(1))
#>       mu_mu    mu_sigma sigma_sigma 
#> -12.9931856   0.9716737  -1.2084218 

# The mean proportion and the dispersion are not orthogonal: the mixed
# entry does not vanish, so their estimates are asymptotically correlated.
eh$mu_sigma[1]
#> [1] 0.9716737

# The strategy argument is inert, the expectation being an exact sum.
identical(eh, distrib_expected_hessian(d, 0:10, th, approx = "mc",
                                       nsim = 50))
#> [1] TRUE
```
