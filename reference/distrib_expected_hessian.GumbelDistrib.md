# Gumbel Expected Hessian

Returns the expectation of the observed Hessian under the model, in
closed form and with no quadrature or simulation. With \\\gamma\\ the
Euler-Mascheroni constant, \$\$\mathbb{E}\left\[\ell^{(\mu\mu)}\right\]
= -\dfrac{1}{\sigma^2}, \qquad
\mathbb{E}\left\[\ell^{(\mu\sigma)}\right\] = \dfrac{1 -
\gamma}{\sigma^2}, \qquad \mathbb{E}\left\[\ell^{(\sigma\sigma)}\right\]
= -\dfrac{(1-\gamma)^2 + \pi^2/6}{\sigma^2}.\$\$

They are closed form for one reason, and it is the same reason the
Weibull's are: under the model \\w = e^{-Z}\\ is standard exponential
whatever the parameters, so every expectation the family needs is a
derivative of \\\Gamma\\ at 2. Here that is \\\mathbb{E}\[w\] = 1\\,
\\\mathbb{E}\[w\log w\] = 1 - \gamma\\ and \\\mathbb{E}\[w(\log w)^2\] =
(1-\gamma)^2 + \pi^2/6 - 1\\, with \\\mathbb{E}\[Z\] = \gamma\\.

**The mixed entry does not vanish.** The location and the scale are not
orthogonal here, where in a symmetric location-scale family such as the
Gaussian they are, and the reason is that this density is skewed. Their
maximum likelihood estimates are asymptotically correlated.

Because the values do not depend on the data, `approx` and `nsim` are
ignored. `y` is read only for its length.

## Arguments

- distrib:

  A `GumbelDistrib` object, from
  [`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. `mu` is not read. `sigma` must be
  strictly positive.

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

## Value

A named list of three numeric vectors, `mu_mu`, `sigma_sigma` and
`mu_sigma`, in that order, each of length `length(y)` and constant
within itself when the parameters are.

## Notation

The **expected information** is
\\\mathbb{E}\[-\partial^2\ell/\partial\theta\\\partial\theta^\top\]\\,
the expectation of the **observed information** under the model. The
Gumbel is a regular family, so the second Bartlett identity holds and
this equals the variance of the score. \\\gamma\\ is the
Euler-Mascheroni constant, \\-\psi(1) \approx 0.5772\\.

## See also

[`distrib_hessian.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GumbelDistrib.md)
for the observed quantity this is the expectation of,
[`distrib_expected_hessian.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Weibull1Distrib.md)
for the family that shares these expectations,
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which inverts it at each step, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- gumbel_distrib()
th <- list(mu = 0, sigma = 1)
e <- distrib_expected_hessian(d, c(-1, 0, 1), th)
lapply(e, unique)
#> $mu_mu
#> [1] -1
#> 
#> $sigma_sigma
#> [1] -1.823681
#> 
#> $mu_sigma
#> [1] 0.4227843
#> 

# Written out with the Euler-Mascheroni constant.
g <- -digamma(1)
c(-1, -((1 - g)^2 + pi^2 / 6), 1 - g)
#> [1] -1.0000000 -1.8236807  0.4227843

# The mixed entry is not zero, so the location and the scale are correlated.
e$mu_sigma[1]
#> [1] 0.4227843

# Negative definite all the same.
M <- matrix(c(e$mu_mu[1], e$mu_sigma[1], e$mu_sigma[1], e$sigma_sigma[1]), 2)
eigen(M, only.values = TRUE)$values
#> [1] -0.8216208 -2.0020598

# The observed Hessian averages onto them over a large sample.
set.seed(11)
z <- distrib_rng(d, 2e5, th)
rbind(observed = vapply(distrib_hessian(d, z, th), mean, numeric(1)),
      expected = vapply(distrib_expected_hessian(d, z, th),
                        function(v) v[1], numeric(1)))
#>               mu_mu sigma_sigma  mu_sigma
#> observed -0.9976522   -1.813257 0.4166854
#> expected -1.0000000   -1.823681 0.4227843
```
