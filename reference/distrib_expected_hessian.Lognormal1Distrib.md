# Lognormal Expected Hessian

Returns the expectation of the observed Hessian under the model, in
closed form and with no quadrature or simulation:
\$\$\mathbb{E}\left\[\ell^{(\mu\mu)}\right\] = -\dfrac{1}{\sigma^2},
\qquad \mathbb{E}\left\[\ell^{(\sigma^2\sigma^2)}\right\] =
-\dfrac{1}{2\sigma^4}, \qquad
\mathbb{E}\left\[\ell^{(\mu\sigma^2)}\right\] = 0.\$\$ They follow from
\\\mathbb{E}\[\log Y\] = \mu\\ and \\\mathbb{E}\[(\log Y - \mu)^2\] =
\sigma^2\\, \\\log Y\\ being Gaussian by construction. The zero
off-diagonal says the two parameters are orthogonal, so their estimates
are asymptotically independent.

The information does not depend on the mean of \\Y\\ at all, only on the
variance of its logarithm. Because the values do not depend on the data,
`approx` and `nsim` are ignored and `y` is read only for its length.

## Arguments

- distrib:

  A `Lognormal1Distrib` object, from
  [`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with components `mu` and `sigma2`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma2` must be strictly positive.

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

A named list of three numeric vectors, `mu_mu`, `sigma2_sigma2` and
`mu_sigma2`, in that order, each of length
`max(length(y), length(mu), length(sigma2))` and constant within itself
when the parameters are.

## Notation

The **expected information** is
\\\mathbb{E}\[-\partial^2\ell/\partial\theta\\\partial\theta^\top\]\\,
the expectation of the **observed information** under the model. The
lognormal is a regular family, so the second Bartlett identity holds and
this equals the variance of the score.

## See also

[`distrib_hessian.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Lognormal1Distrib.md)
for the observed quantity this is the expectation of,
[`distrib_expected_hessian.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gaussian2Distrib.md),
which returns the same numbers,
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which inverts it at each step, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- lognormal1_distrib()
th <- list(mu = 0.5, sigma2 = 0.36)

# The three constants, one value per observation.
lapply(distrib_expected_hessian(d, c(0.5, 1.6, 4), th), unique)
#> $mu_mu
#> [1] -2.777778
#> 
#> $sigma2_sigma2
#> [1] -3.858025
#> 
#> $mu_sigma2
#> [1] 0
#> 
c(-1 / 0.36, -1 / (2 * 0.36^2), 0)
#> [1] -2.777778 -3.858025  0.000000

# The observed Hessian averages onto them over a large sample.
set.seed(11)
z <- distrib_rng(d, 2e4, th)
rbind(observed = vapply(distrib_hessian(d, z, th), mean, numeric(1)),
      expected = vapply(distrib_expected_hessian(d, z, th),
                        function(v) v[1], numeric(1)))
#>              mu_mu sigma2_sigma2   mu_sigma2
#> observed -2.777778     -3.808307 -0.01203833
#> expected -2.777778     -3.858025  0.00000000

# The information does not move with the mean of Y, only with sigma2.
vapply(c(-2, 0, 5),
       function(m) distrib_expected_hessian(d, 0,
                     list(mu = m, sigma2 = 0.36))$sigma2_sigma2,
       numeric(1))
#> [1] -3.858025 -3.858025 -3.858025
```
