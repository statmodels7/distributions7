# Laplace Expected Hessian

Returns the negative of the Fisher information, in closed form and with
no quadrature or simulation: \$\$-I(\mu) = -\dfrac{1}{\sigma^2}, \qquad
-I(\sigma) = -\dfrac{1}{\sigma^2}, \qquad -I(\mu, \sigma) = 0.\$\$

For \\\mu\\ this is **not** the expectation of
[`distrib_hessian.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.LaplaceDistrib.md),
which is identically zero. It is the variance of the score:
\\\partial\ell/\partial\mu = \mathrm{sign}(r)/\sigma\\ has mean 0 and
square \\1/\sigma^2\\ almost surely, so \\\mathrm{Var} = 1/\sigma^2\\.
The mixed entry vanishes because the family is symmetric about \\\mu\\,
so \\\mathbb{E}\[\mathrm{sign}(r)\] = 0\\.

Because the values do not depend on the data, `approx` and `nsim` are
ignored. `y` is read only for its length.

## Arguments

- distrib:

  A `LaplaceDistrib` object, from
  [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored, a closed form being available. Accepted so that the signature
  matches the generic's. This matters more here than elsewhere: the
  strategies that average the observed second derivative would return 0
  for the location, and only the score-based one recovers
  \\1/\sigma^2\\.

- nsim:

  Ignored, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of three numeric vectors, `mu_mu`, `sigma_sigma` and
`mu_sigma`, each of length `length(y)` and constant within itself when
the parameters are.

## Notation

The **observed information** is
\\-\partial^2\ell/\partial\theta\\\partial\theta^\top\\ at the data. The
**expected information** is usually its expectation, and for a regular
family the two agree with the variance of the score by the second
Bartlett identity. The Laplace is **not** regular in \\\mu\\: the first
identity \\\mathbb{E}\[\partial\ell/\partial\mu\] = 0\\ still holds, the
second does not, and the information is **defined** as the variance of
the score. That is what this method returns.

## See also

[`distrib_hessian.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.LaplaceDistrib.md),
whose `mu_mu` is zero;
[`distrib_gradient.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.LaplaceDistrib.md),
whose variance this is;
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which inverts this matrix and so fits this family where a Newton step
could not; and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic and the strategies it offers a family with no closed
form.

## Examples

``` r
d <- laplace_distrib()
th <- list(mu = 0.4, sigma = 1.5)

# Both diagonal entries are -1/sigma^2; the mixed entry is 0.
lapply(distrib_expected_hessian(d, c(-1.2, 0.3, 2.5), th), unique)
#> $mu_mu
#> [1] -0.4444444
#> 
#> $sigma_sigma
#> [1] -0.4444444
#> 
#> $mu_sigma
#> [1] 0
#> 
-1 / 1.5^2
#> [1] -0.4444444

# It is the variance of the score, which the sample confirms; the mean of
# the observed second derivative is 0 and does not.
set.seed(12)
z <- distrib_rng(d, 1e5, th)
s <- distrib_gradient(d, z, th)$mu
c(var_of_score = mean(s^2),
  information = -distrib_expected_hessian(d, 0, th)$mu_mu,
  mean_observed = mean(distrib_hessian(d, z, th)$mu_mu))
#>  var_of_score   information mean_observed 
#>     0.4444444     0.4444444     0.0000000 

# Fisher scoring can fit the family because this matrix is nonsingular;
# a Newton step on the observed Hessian would divide by zero in mu.
coef(fit_distrib(d, z))
#>        mu     sigma 
#> 0.4051025 1.4992390 
```
