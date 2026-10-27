# Gaussian Expected Hessian

Returns the expectation of the observed Hessian under the model, in
closed form and with no quadrature or simulation:
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{1}{\sigma^2}, \qquad \mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \sigma^2}\right\] = -\dfrac{2}{\sigma^2}, \qquad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu \\ \partial
\sigma}\right\] = 0.\$\$ They follow from \\\mathbb{E}\[(Y-\mu)^2\] =
\sigma^2\\ and \\\mathbb{E}\[Y-\mu\] = 0\\. The negative of this matrix
is the Fisher information for one observation. Its off-diagonal entry is
zero, so the mean and the standard deviation are orthogonal and their
estimates are asymptotically independent.

Because the values do not depend on the data, `y` is read only for its
length.

## Arguments

- distrib:

  A `Gaussian1Distrib` object, from
  [`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` must be strictly positive.

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

A named list of three numeric vectors, `mu_mu`, `sigma_sigma` and
`mu_sigma`, each of length `max(length(y), length(mu), length(sigma))`
and constant within itself when the parameters are.

## Notation

The expected information is
\\\mathbb{E}\[-\partial^2\ell/\partial\theta\\\partial\theta^\top\]\\,
the expectation of the observed information under the model. The
Gaussian is a regular family, so the second Bartlett identity holds and
this equals the variance of the score.

## See also

[`distrib_hessian.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Gaussian1Distrib.md)
for the observed quantity whose expectation this is,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
and
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which invert it at each step, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- gaussian1_distrib()
th <- list(mu = 0.4, sigma = 1.5)

# The three components are constants, with one value per observation.
lapply(distrib_expected_hessian(d, c(-1.2, 0.3, 2.5), th), unique)
#> $mu_mu
#> [1] -0.4444444
#> 
#> $sigma_sigma
#> [1] -0.8888889
#> 
#> $mu_sigma
#> [1] 0
#> 

# The observed Hessian averages onto them over a large sample.
set.seed(11)
z <- distrib_rng(d, 2e4, th)
rbind(observed = vapply(distrib_hessian(d, z, th), mean, numeric(1)),
      expected = vapply(distrib_expected_hessian(d, z, th),
                        function(v) v[1], numeric(1)))
#>               mu_mu sigma_sigma     mu_sigma
#> observed -0.4444444  -0.8802976 -0.002311359
#> expected -0.4444444  -0.8888889  0.000000000

# The mean and the standard deviation are orthogonal: the mixed entry is 0.
distrib_expected_hessian(d, 0, th)$mu_sigma
#> [1] 0
```
