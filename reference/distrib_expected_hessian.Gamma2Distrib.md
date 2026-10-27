# Gamma Expected Hessian in Mean and Variance

Returns the expectation of the observed Hessian under the model, in
closed form and with no quadrature or simulation. With \\\alpha =
\mu^2/\sigma^2\\ and \\\psi_1\\ the trigamma function,
\$\$\mathbb{E}\left\[\ell^{(\mu\mu)}\right\] = \dfrac{3\sigma^2 -
4\mu^2\psi_1(\alpha)}{(\sigma^2)^2}, \qquad
\mathbb{E}\left\[\ell^{(\sigma^2\sigma^2)}\right\] =
-\dfrac{\mu^2\left\\\mu^2\psi_1(\alpha) - \sigma^2\right\\}
{(\sigma^2)^4},\$\$ \$\$\mathbb{E}\left\[\ell^{(\mu\sigma^2)}\right\] =
\dfrac{2\mu\left\\\mu^2\psi_1(\alpha) - \sigma^2\right\\}
{(\sigma^2)^3}.\$\$ They follow from \\\mathbb{E}\[Y\] = \mu\\ and
\\\mathbb{E}\[\log Y\] = \psi(\alpha) - \log\lambda\\.

The mixed entry does not vanish, so the mean and the variance are not
orthogonal in this parametrization; their maximum likelihood estimates
are asymptotically correlated.
[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md)
carries the same law with a dispersion in place of the variance, and
there the mixed entry is exactly zero. The matrix is negative definite
throughout, so the information is positive definite.

The kernel evaluates them through the decomposition described on
[`distrib_deriv3.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Gamma2Distrib.md),
which keeps the digits at a large shape.

Because the values do not depend on the data, `y` is read only for its
length.

## Arguments

- distrib:

  A `Gamma2Distrib` object, from
  [`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with components `mu` and `sigma2`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). The
  generic reads it, and this method does not.

- method:

  `NULL` or an object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md).
  The generic applies an object before dispatching on the family, so
  this method receives `NULL`, and it returns the closed form of the
  family. See
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  a fixed number of observations, set per kernel, the kernel runs
  sequentially. Defaults to `1L`.

## Value

A named list of three numeric vectors, `mu_mu`, `sigma2_sigma2` and
`mu_sigma2`, in that order, each of length
`max(length(y), length(mu), length(sigma2))` and constant within itself
when the parameters are.

## Notation

The expected information is
\\\mathbb{E}\[-\partial^2\ell/\partial\theta\\\partial\theta^\top\]\\,
the expectation of the observed information under the model. The gamma
is a regular family, so the second Bartlett identity holds and this
equals the variance of the score.

## See also

[`distrib_hessian.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Gamma2Distrib.md)
for the observed quantity whose expectation this is,
[`distrib_expected_hessian.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gamma1Distrib.md)
for the orthogonal parametrization of the same law,
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which inverts it at each step, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- gamma2_distrib()
th <- list(mu = 3, sigma2 = 2)

# The three entries are constants, with one value per observation.
eh <- lapply(distrib_expected_hessian(d, c(1, 3, 5), th), unique)
eh
#> $mu_mu
#> [1] -0.7385259
#> 
#> $sigma2_sigma2
#> [1] -0.1341708
#> 
#> $mu_sigma2
#> [1] 0.1788944
#> 

# The three formulas are written out with the trigamma function.
al <- 9 / 2
c((3 * 2 - 4 * 9 * trigamma(al)) / 2^2,
  -9 * (9 * trigamma(al) - 2) / 2^4,
  2 * 3 * (9 * trigamma(al) - 2) / 2^3)
#> [1] -0.7385259 -0.1341708  0.1788944

# The mixed entry is not zero: the mean and the variance are correlated.
# In gamma1, the same law in the dispersion, it is exactly zero.
c(gamma2 = eh$mu_sigma2,
  gamma1 = distrib_expected_hessian(gamma1_distrib(), 0,
                                    list(mu = 3, phi = 2 / 9))$mu_phi)
#>    gamma2    gamma1 
#> 0.1788944 0.0000000 

# Negative definite all the same, so the information is positive definite.
M <- matrix(c(eh$mu_mu, eh$mu_sigma2, eh$mu_sigma2, eh$sigma2_sigma2), 2)
eigen(M, only.values = TRUE)$values
#> [1] -0.08518675 -0.78751001

# The observed Hessian averages onto them over a large sample.
set.seed(11)
z <- distrib_rng(d, 2e4, th)
rbind(observed = vapply(distrib_hessian(d, z, th), mean, numeric(1)),
      expected = vapply(distrib_expected_hessian(d, z, th),
                        function(v) v[1], numeric(1)))
#>               mu_mu sigma2_sigma2 mu_sigma2
#> observed -0.7417375    -0.1381203 0.1826196
#> expected -0.7385259    -0.1341708 0.1788944
```
