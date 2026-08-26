# Cauchy Observed Hessian

Computes the three distinct second derivatives of the Cauchy log-density
with respect to \\\mu\\ and \\\sigma\\, one value per observation, in
closed form. Writing \\r = y - \mu\\ and \\d = \sigma^2 + r^2\\,
\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{2(r^2 -
\sigma^2)}{d^2}, \qquad \dfrac{\partial^2 \ell}{\partial \sigma^2} =
\dfrac{\sigma^4 - 4\sigma^2 r^2 - r^4}{\sigma^2 d^2}, \qquad
\dfrac{\partial^2 \ell}{\partial \mu \\ \partial \sigma} =
-\dfrac{4\sigma r}{d^2}.\$\$

The curvature in \\\mu\\ turns **positive** wherever \\\|r\| \>
\sigma\\, so a single observation beyond one scale unit contributes
convexity and the observed information can fail to be positive definite.
The expected values are well behaved, which is why Fisher scoring is the
steadier route on this family; see
[`distrib_expected_hessian.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.CauchyDistrib.md).

## Arguments

- distrib:

  A `CauchyDistrib` object, from
  [`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count says. Defaults to `1L`.

## Value

A named list of three numeric vectors, `mu_mu`, `sigma_sigma` and
`mu_sigma`, each of length `max(length(y), length(mu), length(sigma))`.
The three name the distinct entries of a symmetric \\2 \times 2\\ matrix
per observation.

## See also

[`distrib_gradient.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.CauchyDistrib.md)
for the score,
[`distrib_expected_hessian.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.CauchyDistrib.md)
for the expectation of this quantity,
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md)
for the estimation route that uses it, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- cauchy_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)
h <- distrib_hessian(d, y, th)

# The closed forms, written out.
r <- y - 0.4; dd <- 1.5^2 + r^2
all.equal(h$mu_mu, 2 * (r^2 - 1.5^2) / dd^2)
#> [1] TRUE
all.equal(h$mu_sigma, -4 * 1.5 * r / dd^2)
#> [1] TRUE

# The curvature in mu is positive wherever |r| exceeds sigma.
data.frame(r = r, mu_mu = h$mu_mu, beyond_sigma = abs(r) > 1.5)
#>      r       mu_mu beyond_sigma
#> 1 -1.6  0.02679795         TRUE
#> 2 -0.1 -0.87712429        FALSE
#> 3  2.1  0.09739469         TRUE

# A central difference of the score reproduces it.
eps <- 1e-5
up <- distrib_gradient(d, y, list(mu = 0.4 + eps, sigma = 1.5))$mu
dn <- distrib_gradient(d, y, list(mu = 0.4 - eps, sigma = 1.5))$mu
all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-6)
#> [1] TRUE
```
