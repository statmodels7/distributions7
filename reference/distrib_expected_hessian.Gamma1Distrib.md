# Gamma Expected Hessian in Mean and Dispersion

Returns the expectation of the observed Hessian under the model, in
closed form and with no quadrature or simulation. With \\s = 1/\phi\\,
\$\$\mathbb{E}\left\[\ell^{(\mu\mu)}\right\] = -\dfrac{1}{\phi\mu^2},
\qquad \mathbb{E}\left\[\ell^{(\mu\phi)}\right\] = 0, \qquad
\mathbb{E}\left\[\ell^{(\phi\phi)}\right\] = s^4\left\\\dfrac{1}{s} -
\psi'(s)\right\\.\$\$ They follow from \\\mathbb{E}\[Y\] = \mu\\ and
\\\mathbb{E}\[\log(Y/\mu)\] = \psi(s) - \log s\\, the second of which is
exactly what makes the score in \\\phi\\ have mean zero. The quantity
\\1/s - \psi'(s)\\ is negative for every \\s \> 0\\, so the pure
dispersion entry is negative and the information is positive definite
everywhere, where the observed Hessian is not.

The zero off-diagonal says the mean and the dispersion are orthogonal.
That orthogonality is the reason this parametrization is the one a
generalized linear model uses: the mean equation can be fitted with the
dispersion held at any value without biasing it.

Because the values do not depend on the data, `approx` and `nsim` are
ignored. `y` is read only for its length.

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

A named list of three numeric vectors, `mu_mu`, `mu_phi` and `phi_phi`,
each of length `max(length(y), length(mu), length(phi))` and constant
within itself when the parameters are.

## Notation

The **expected information** is
\\\mathbb{E}\[-\partial^2\ell/\partial\theta\\\partial\theta^\top\]\\,
the expectation of the **observed information** under the model. The
gamma is a regular family, so the second Bartlett identity holds and
this equals the variance of the score. \\\psi'\\ is the trigamma
function.

## See also

[`distrib_hessian.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Gamma1Distrib.md)
for the observed quantity this is the expectation of,
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which inverts it at each step, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- gamma1_distrib()
th <- list(mu = 3, phi = 0.5)

# The three constants, one value per observation.
lapply(distrib_expected_hessian(d, c(1, 3, 5), th), unique)
#> $mu_mu
#> [1] -0.2222222
#> 
#> $mu_phi
#> [1] 0
#> 
#> $phi_phi
#> [1] -2.318945
#> 

# The dispersion entry, written out with the trigamma function.
s <- 1 / 0.5
s^4 * (1 / s - trigamma(s))
#> [1] -2.318945

# The observed Hessian averages onto them over a large sample.
set.seed(11)
z <- distrib_rng(d, 2e4, th)
rbind(observed = vapply(distrib_hessian(d, z, th), mean, numeric(1)),
      expected = vapply(distrib_expected_hessian(d, z, th),
                        function(v) v[1], numeric(1)))
#>               mu_mu        mu_phi   phi_phi
#> observed -0.2222366 -4.310552e-05 -2.288278
#> expected -0.2222222  0.000000e+00 -2.318945

# Negative definite at every parameter setting, where the observed Hessian
# is positive in mu wherever y < mu/2.
phi <- c(0.05, 0.5, 5)
vapply(phi, function(p) {
  s <- 1 / p
  s^4 * (1 / s - trigamma(s))
}, numeric(1))
#> [1] -203.3316696   -2.3189451   -0.0340278
```
