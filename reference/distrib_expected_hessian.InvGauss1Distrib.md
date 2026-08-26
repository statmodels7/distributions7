# Inverse Gaussian Expected Hessian in Mean and Dispersion

Returns the expectation of the observed Hessian under the model, in
closed form and with no quadrature or simulation:
\$\$\mathbb{E}\left\[\ell^{(\mu\mu)}\right\] = -\dfrac{1}{\phi\mu^3},
\qquad \mathbb{E}\left\[\ell^{(\phi\phi)}\right\] = -\dfrac{1}{2\phi^2},
\qquad \mathbb{E}\left\[\ell^{(\mu\phi)}\right\] = 0.\$\$ They follow
from \\\mathbb{E}\[Y\] = \mu\\ and \\\mathbb{E}\[(Y-\mu)^2/Y\] =
\phi\mu^2\\. The second identity is the one that gives the score in
\\\phi\\ mean zero.

Both diagonal entries are negative at every parameter setting, so the
information is positive definite everywhere, where the observed Hessian
of
[`distrib_hessian.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.InvGauss1Distrib.md)
is not. The zero off-diagonal says the mean and the dispersion are
orthogonal, so the mean equation can be fitted with the dispersion held
at any value without biasing it: this is the generalized linear model
parametrization, with variance function \\V(\mu) = \mu^3\\.

Because the values do not depend on the data, `approx` and `nsim` are
ignored. `y` is read only for its length.

## Arguments

- distrib:

  An `InvGauss1Distrib` object, from
  [`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with components `mu` and `phi`, each a numeric vector of
  length 1 or of the length of `y`. A component of length 1 is recycled.
  Both must be strictly positive.

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
inverse Gaussian is a regular family, so the second Bartlett identity
holds and this equals the variance of the score.

## See also

[`distrib_hessian.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.InvGauss1Distrib.md)
for the observed quantity this is the expectation of,
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which inverts it at each step, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- invgauss1_distrib()
th <- list(mu = 1, phi = 2)

# The three constants, one value per observation.
lapply(distrib_expected_hessian(d, c(0.5, 1, 2), th), unique)
#> $mu_mu
#> [1] -0.5
#> 
#> $phi_phi
#> [1] -0.125
#> 
#> $mu_phi
#> [1] 0
#> 
c(-1 / (2 * 1^3), -1 / (2 * 2^2), 0)
#> [1] -0.500 -0.125  0.000

# Negative definite at every setting, where the observed Hessian is
# positive in phi at y = mu.
c(expected = distrib_expected_hessian(d, 1, th)$phi_phi,
  observed = distrib_hessian(d, 1, th)$phi_phi)
#> expected observed 
#>   -0.125    0.125 

# The observed Hessian averages onto them over a large sample.
set.seed(11)
z <- distrib_rng(d, 2e4, th)
rbind(observed = vapply(distrib_hessian(d, z, th), mean, numeric(1)),
      expected = vapply(distrib_expected_hessian(d, z, th),
                        function(v) v[1], numeric(1)))
#>               mu_mu    phi_phi      mu_phi
#> observed -0.5235085 -0.1233891 -0.00391808
#> expected -0.5000000 -0.1250000  0.00000000
```
