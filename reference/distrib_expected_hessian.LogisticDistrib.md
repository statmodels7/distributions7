# Logistic Expected Hessian

Returns the expectation of the observed Hessian under the model, in
closed form and with no quadrature or simulation:
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{1}{3\sigma^2}, \qquad \mathbb{E}\left\[\dfrac{\partial^2
\ell}{\partial \sigma^2}\right\] = -\dfrac{3 + \pi^2}{9\sigma^2}, \qquad
\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu \\ \partial
\sigma}\right\] = 0.\$\$ The mixed entry vanishes because the family is
symmetric about \\\mu\\, so the location and the scale are orthogonal
and their estimates are asymptotically independent. The information in
the location is \\1/(3\sigma^2)\\, a third of a Gaussian's
\\1/\sigma^2\\ at the same \\\sigma\\. Compared at equal variance the
two are much closer, the logistic variance being \\\pi^2\sigma^2/3\\.

Because the values do not depend on the data, `approx` and `nsim` are
ignored. `y` is read only for its length.

## Arguments

- distrib:

  A `LogisticDistrib` object, from
  [`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md).

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
logistic is a regular family, so the second Bartlett identity holds and
this equals the variance of the score.

## See also

[`distrib_hessian.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.LogisticDistrib.md)
for the observed quantity this is the expectation of,
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which inverts it at each step, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- logistic_distrib()
th <- list(mu = 0.4, sigma = 1.5)

# The three closed forms.
lapply(distrib_expected_hessian(d, c(-1.2, 0.3, 2.5), th), unique)
#> $mu_mu
#> [1] -0.1481481
#> 
#> $sigma_sigma
#> [1] -0.635536
#> 
#> $mu_sigma
#> [1] 0
#> 
c(-1 / (3 * 1.5^2), -(3 + pi^2) / (9 * 1.5^2))
#> [1] -0.1481481 -0.6355360

# The observed Hessian averages onto them over a large sample.
set.seed(8)
z <- distrib_rng(d, 2e5, th)
rbind(observed = vapply(distrib_hessian(d, z, th), mean, numeric(1)),
      expected = vapply(distrib_expected_hessian(d, z, th),
                        function(v) v[1], numeric(1)))
#>               mu_mu sigma_sigma      mu_sigma
#> observed -0.1480978  -0.6370314 -0.0003646557
#> expected -0.1481481  -0.6355360  0.0000000000
```
