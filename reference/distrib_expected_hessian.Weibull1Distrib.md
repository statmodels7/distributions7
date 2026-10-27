# Weibull Expected Hessian

Returns the expectation of the observed Hessian under the model, in
closed form and with no quadrature or simulation:
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{\sigma^2}{\mu^2}, \qquad \mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \sigma^2}\right\] =
-\dfrac{1}{\sigma^2}\left\\(1-\gamma)^2 + \dfrac{\pi^2}{6}\right\\,
\qquad \mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu \\ \partial
\sigma}\right\] = \dfrac{1 - \gamma}{\mu},\$\$ with \\\gamma \method
0.5772\\ the Euler-Mascheroni constant. The mixed entry does not vanish,
so the scale and the shape are not orthogonal and their estimates are
asymptotically correlated.

Because a closed form exists, the method returns the same three numbers
whatever the call.

## Arguments

- distrib:

  A `Weibull1Distrib` object, from
  [`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

- y:

  A numeric vector of observations. Only its length is read, because the
  expectation does not depend on the data; the values are ignored.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
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

A named list of three numeric vectors, `mu_mu`, `sigma_sigma` and
`mu_sigma`, each of length `max(length(y), length(mu), length(sigma))`
and constant when `mu` and `sigma` have length 1.

## Why the expectations are gamma derivatives

Under the model \\u = (Y/\mu)^{\sigma}\\ is standard exponential
whatever the parameters are, so every expectation that the Hessian needs
is a moment of \\u\\ against a power of \\\log u\\, and each of those is
a derivative of \\\Gamma\\ at 2: \$\$\mathbb{E}\[u\] = 1, \qquad
\mathbb{E}\[u \log u\] = 1 - \gamma, \qquad \mathbb{E}\[u (\log u)^2\] =
(1-\gamma)^2 + \dfrac{\pi^2}{6} - 1.\$\$ Substituting \\\log z = (\log
u)/\sigma\\ into the three observed components gives the display above.
The Gumbel family shares the substitution, since
\\\exp(-\text{Gumbel})\\ is Weibull, so
[`distrib_expected_hessian.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GumbelDistrib.md)
rests on the same three moments.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu \> 0\\ the scale,
\\\sigma \> 0\\ the shape and \\\gamma\\ the Euler-Mascheroni constant,
`-digamma(1)`.

## See also

[`distrib_hessian.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Weibull1Distrib.md)
for the observed Hessian, of which this is the expectation,
[`distrib_gradient.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Weibull1Distrib.md)
for the score, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- weibull1_distrib()
y <- c(0.5, 1.2, 3.0)
th <- list(mu = 2, sigma = 1.5)
eh <- distrib_expected_hessian(d, y, th)
vapply(eh, function(z) z[1], numeric(1))
#>       mu_mu sigma_sigma    mu_sigma 
#>  -0.5625000  -0.8105247   0.2113922 

# The closed forms are written out.
eg <- -digamma(1)
c(mu_mu = -1.5^2 / 4,
  sigma_sigma = -((1 - eg)^2 + pi^2 / 6) / 1.5^2,
  mu_sigma = (1 - eg) / 2)
#>       mu_mu sigma_sigma    mu_sigma 
#>  -0.5625000  -0.8105247   0.2113922 

# Averaging the observed Hessian over draws reaches the same three numbers.
set.seed(4)
z <- distrib_rng(d, 2e5, th)
vapply(distrib_hessian(d, z, th), mean, numeric(1))
#>       mu_mu sigma_sigma    mu_sigma 
#>  -0.5635714  -0.8109912   0.2125961 

# The method argument has no effect, because the expectation is exact.
identical(eh, distrib_expected_hessian(d, y, th, method = expect_mc(nsim = 50)))
#> [1] FALSE
```
