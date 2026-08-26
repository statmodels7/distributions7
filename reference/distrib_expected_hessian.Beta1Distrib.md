# Beta Expected Hessian in Mean and Precision

Returns the expectation of the observed Hessian under the model, in
closed form and with no quadrature or simulation. With \\\alpha =
\mu\phi\\, \\\beta = (1-\mu)\phi\\ and \\\psi_1\\ the trigamma function,
\$\$\mathbb{E}\left\[\ell^{(\mu\mu)}\right\] =
-\phi^2\left\\\psi_1(\alpha) + \psi_1(\beta)\right\\, \qquad
\mathbb{E}\left\[\ell^{(\phi\phi)}\right\] = \psi_1(\phi) -
\mu^2\psi_1(\alpha) - (1-\mu)^2\psi_1(\beta),\$\$
\$\$\mathbb{E}\left\[\ell^{(\mu\phi)}\right\] =
-\phi\left\\\mu\psi_1(\alpha) - (1-\mu)\psi_1(\beta)\right\\.\$\$

The first two are the observed values themselves, being free of the
data. The third drops the log-odds residual the observed mixed entry
carries, whose expectation is zero because
\\\mathbb{E}\[\log\\Y/(1-Y)\\\] = \psi(\alpha) - \psi(\beta)\\.

The mixed entry does not vanish, so the mean and the precision are not
orthogonal; their maximum likelihood estimates are asymptotically
correlated. It does vanish at \\\mu = 1/2\\, where the two shapes are
equal and the density is symmetric.

Because the values do not depend on the data, `approx` and `nsim` are
ignored. `y` is read only for its length.

## Arguments

- distrib:

  A `Beta1Distrib` object, from
  [`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with components `mu` and `phi`, each a numeric vector of
  length 1 or of the length of `y`. A component of length 1 is recycled.
  `mu` must lie strictly in \\(0, 1)\\ and `phi` must be strictly
  positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored, a closed form being available. Accepted so that the signature
  matches the generic's, where it selects between the Bartlett,
  quadrature, Monte Carlo and outer-product routes.

- nsim:

  Ignored, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count says. Defaults to `1L`.

## Value

A named list of three numeric vectors, `mu_mu`, `phi_phi` and `mu_phi`,
in that order, each of length `max(length(y), length(mu), length(phi))`
and constant within itself when the parameters are.

## Notation

The **expected information** is
\\\mathbb{E}\[-\partial^2\ell/\partial\theta\\\partial\theta^\top\]\\,
the expectation of the **observed information** under the model. The
beta is a regular family, so the second Bartlett identity holds and this
equals the variance of the score.

## See also

[`distrib_hessian.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Beta1Distrib.md)
for the observed quantity this is the expectation of,
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which inverts it at each step, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- beta1_distrib()
th <- list(mu = 0.4, phi = 5)
eh <- distrib_expected_hessian(d, c(0.2, 0.5, 0.8), th)
lapply(eh, unique)
#> $mu_mu
#> [1] -25.9967
#> 
#> $phi_phi
#> [1] -0.02404276
#> 
#> $mu_phi
#> [1] -0.1050659
#> 

# Written out with the trigamma function.
a <- 0.4 * 5
b <- 0.6 * 5
c(-5^2 * (trigamma(a) + trigamma(b)),
  trigamma(5) - 0.4^2 * trigamma(a) - 0.6^2 * trigamma(b),
  -5 * (0.4 * trigamma(a) - 0.6 * trigamma(b)))
#> [1] -25.99670334  -0.02404276  -0.10506593

# The observed mixed entry averages onto it over a large sample; the other
# two are equal to it observation by observation.
set.seed(9)
z <- distrib_rng(d, 5e5, th)
c(observed = mean(distrib_hessian(d, z, th)$mu_phi),
  expected = eh$mu_phi[1])
#>   observed   expected 
#> -0.1059716 -0.1050659 

# The mixed entry vanishes at mu = 1/2, where the density is symmetric.
distrib_expected_hessian(d, 0.5, list(mu = 0.5, phi = 5))$mu_phi
#> [1] 0
```
