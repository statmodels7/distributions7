# Gamma Expected Hessian in Mean and Dispersion

Returns the expectation of the observed Hessian under the model, in
closed form and with no quadrature or simulation. With \\s = 1/\phi\\,
\$\$\mathbb{E}\left\[\ell^{(\mu\mu)}\right\] = -\dfrac{1}{\phi\mu^2},
\qquad \mathbb{E}\left\[\ell^{(\mu\phi)}\right\] = 0, \qquad
\mathbb{E}\left\[\ell^{(\phi\phi)}\right\] = s^4\left\\\dfrac{1}{s} -
\psi'(s)\right\\.\$\$ They follow from \\\mathbb{E}\[Y\] = \mu\\ and
\\\mathbb{E}\[\log(Y/\mu)\] = \psi(s) - \log s\\, and the second of
these identities is exactly the condition under which the score in
\\\phi\\ has mean zero. The quantity \\1/s - \psi'(s)\\ is negative for
every \\s \> 0\\, so the pure dispersion entry is negative and the
information is positive definite everywhere, where the observed Hessian
is not.

The off-diagonal entry is zero, so the mean and the dispersion are
orthogonal. The score in \\\mu\\ has the dispersion as a common factor,
so the mean equation can be solved at any fixed value of the dispersion.

Because the values do not depend on the data, `y` is read only for its
length.

## Arguments

- distrib:

  A `Gamma1Distrib` object, from
  [`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with components `mu` and `phi`, each a numeric vector of
  length 1 or of the length of `y`. A component of length 1 is recycled.
  Both must be strictly positive.

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

A named list of three numeric vectors, `mu_mu`, `mu_phi` and `phi_phi`,
each of length `max(length(y), length(mu), length(phi))` and constant
within itself when the parameters are.

## Notation

The expected information is
\\\mathbb{E}\[-\partial^2\ell/\partial\theta\\\partial\theta^\top\]\\,
the expectation of the observed information under the model. The gamma
is a regular family, so the second Bartlett identity holds and this
equals the variance of the score. \\\psi'\\ is the trigamma function.

## See also

[`distrib_hessian.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Gamma1Distrib.md)
for the observed quantity whose expectation this is,
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which inverts it at each step, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- gamma1_distrib()
th <- list(mu = 3, phi = 0.5)

# The three constants take one value per observation.
lapply(distrib_expected_hessian(d, c(1, 3, 5), th), unique)
#> $mu_mu
#> [1] -0.2222222
#> 
#> $phi_phi
#> [1] -2.318945
#> 
#> $mu_phi
#> [1] 0
#> 

# The dispersion entry is written out with the trigamma function.
s <- 1 / 0.5
s^4 * (1 / s - trigamma(s))
#> [1] -2.318945

# The observed Hessian averages onto them over a large sample.
set.seed(11)
z <- distrib_rng(d, 2e4, th)
rbind(observed = vapply(distrib_hessian(d, z, th), mean, numeric(1)),
      expected = vapply(distrib_expected_hessian(d, z, th),
                        function(v) v[1], numeric(1)))
#>               mu_mu   phi_phi        mu_phi
#> observed -0.2222366 -2.288278 -4.310552e-05
#> expected -0.2222222 -2.318945  0.000000e+00

# The expected Hessian is negative definite at every parameter setting,
# whereas the observed Hessian is positive in mu wherever y < mu/2.
phi <- c(0.05, 0.5, 5)
vapply(phi, function(p) {
  s <- 1 / p
  s^4 * (1 / s - trigamma(s))
}, numeric(1))
#> [1] -203.3316696   -2.3189451   -0.0340278
```
