# Cauchy Expected Hessian

Returns the expectation of the observed Hessian under the model, in
closed form and with no quadrature or simulation:
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right\] =
-\dfrac{1}{2\sigma^2}, \qquad \mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \mu \\ \partial \sigma}\right\] = 0.\$\$ The information
in the location is \\1/(2\sigma^2)\\, half of a Gaussian's
\\1/\sigma^2\\ at the same scale, which is the price the heavy tails
charge. The zero off-diagonal says the location and the scale are
orthogonal. The expectations exist even though no moment of the family
does: they are expectations of bounded functions of \\y\\, not of \\y\\.

Because the values do not depend on the data, `approx` and `nsim` are
ignored. `y` is read only for its length.

## Arguments

- distrib:

  A `CauchyDistrib` object, from
  [`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored, a closed form being available. Accepted so that the signature
  matches the generic's.

- nsim:

  Ignored, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count says. Defaults to `1L`.

## Value

A named list of three numeric vectors, `mu_mu`, `sigma_sigma` and
`mu_sigma`, each of length `max(length(y), length(mu), length(sigma))`.

## Notation

The **expected information** is
\\\mathbb{E}\[-\partial^2\ell/\partial\theta\\\partial\theta^\top\]\\,
the expectation of the **observed information** under the model. The
Cauchy is a regular family in both parameters, so the second Bartlett
identity holds and this equals the variance of the score.

## See also

[`distrib_hessian.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.CauchyDistrib.md)
for the observed quantity this is the expectation of,
[`distrib_expected_hessian.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gaussian1Distrib.md)
for the light-tailed comparison,
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which inverts it at each step, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- cauchy_distrib()
th <- list(mu = 0.4, sigma = 1.5)

# Both diagonal entries are -1/(2 sigma^2); the mixed entry is 0.
lapply(distrib_expected_hessian(d, c(-1.2, 0.3, 2.5), th), unique)
#> $mu_mu
#> [1] -0.2222222
#> 
#> $sigma_sigma
#> [1] -0.2222222
#> 
#> $mu_sigma
#> [1] 0
#> 
-1 / (2 * 1.5^2)
#> [1] -0.2222222

# Half the information a Gaussian of the same scale carries in its location.
distrib_expected_hessian(gaussian1_distrib(), 0, th)$mu_mu
#> [1] -0.4444444

# The observed Hessian averages onto it over a large sample, even though
# the sample mean of the same draws does not converge at all.
set.seed(3)
z <- distrib_rng(d, 2e5, th)
vapply(distrib_hessian(d, z, th), mean, numeric(1))
#>         mu_mu   sigma_sigma      mu_sigma 
#> -0.2225084173 -0.2219360272  0.0001745913 
```
