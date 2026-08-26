# Gaussian Expected Hessian in Mean and Precision

Returns the expectation of the observed Hessian under the model, in
closed form and with no quadrature or simulation:
\$\$\mathbb{E}\left\[\ell^{(\mu\mu)}\right\] = -\tau, \qquad
\mathbb{E}\left\[\ell^{(\mu\tau)}\right\] = 0, \qquad
\mathbb{E}\left\[\ell^{(\tau\tau)}\right\] = -\dfrac{1}{2\tau^2}.\$\$
They follow from \\\mathbb{E}\[Y-\mu\] = 0\\, the other two entries
being free of the data already. The negative of this matrix is the
Fisher information for one observation; the zero off-diagonal says the
mean and the precision are orthogonal, so their estimates are
asymptotically independent.

Because the values do not depend on the data, `approx` and `nsim` are
ignored. `y` is read only for its length.

## Arguments

- distrib:

  A `Gaussian3Distrib` object, from
  [`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with components `mu` and `tau`, each a numeric vector of
  length 1 or of the length of `y`. A component of length 1 is recycled.
  `tau` must be strictly positive.

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

A named list of three numeric vectors, `mu_mu`, `mu_tau` and `tau_tau`,
each of length `max(length(y), length(mu), length(tau))` and constant
within itself when the parameters are.

## Notation

The **expected information** is
\\\mathbb{E}\[-\partial^2\ell/\partial\theta\\\partial\theta^\top\]\\,
the expectation of the **observed information** under the model. The
Gaussian is a regular family, so the second Bartlett identity holds and
this equals the variance of the score.

## See also

[`distrib_hessian.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Gaussian3Distrib.md)
for the observed quantity this is the expectation of,
[`distrib_expected_hessian.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gaussian1Distrib.md)
for the same information in the standard deviation,
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which inverts it at each step, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- gaussian3_distrib()
th <- list(mu = 1, tau = 0.25)

# The three constants, one value per observation.
lapply(distrib_expected_hessian(d, c(-1.2, 0.3, 2.5), th), unique)
#> $mu_mu
#> [1] -0.25
#> 
#> $mu_tau
#> [1] 0
#> 
#> $tau_tau
#> [1] -8
#> 

# The observed Hessian averages onto them over a large sample.
set.seed(11)
z <- distrib_rng(d, 2e4, th)
rbind(observed = vapply(distrib_hessian(d, z, th), mean, numeric(1)),
      expected = vapply(distrib_expected_hessian(d, z, th),
                        function(v) v[1], numeric(1)))
#>          mu_mu      mu_tau tau_tau
#> observed -0.25 0.005200559      -8
#> expected -0.25 0.000000000      -8

# The information transforms by the delta method: with dtau/dsigma =
# -2/sigma^3, the precision entry carries onto gaussian1's -2/sigma^2.
distrib_expected_hessian(d, 0, th)$tau_tau * (2 / 2^3)^2
#> [1] -0.5
distrib_expected_hessian(gaussian1_distrib(), 0,
                         list(mu = 1, sigma = 2))$sigma_sigma
#> [1] -0.5
```
