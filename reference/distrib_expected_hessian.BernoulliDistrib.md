# Bernoulli Expected Hessian

Returns the expectation of the observed second derivative under the
model, in closed form and with no quadrature or simulation:
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{1}{\mu(1-\mu)},\$\$ which follows from \\\mathbb{E}\[Y\] = \mu\\
in \\-y/\mu^2 - (1-y)/(1-\mu)^2\\. The Fisher information for one
observation is \\1/(\mu(1-\mu))\\, the reciprocal of the variance, so it
is smallest at \\\mu = 1/2\\ and grows without bound as the probability
approaches either endpoint.

On the **link** scale with the default logit the value is
\\-\mu(1-\mu)\\, and so is the observed Hessian: the logit is the
canonical link, the observed curvature there carries no data, and the
two coincide exactly. Fisher scoring and Newton's method therefore take
the same step on a logistic regression, so iteratively reweighted least
squares is both at once.

Because the value does not depend on the data, `approx` and `nsim` are
ignored. `y` is read only for its length.

## Arguments

- distrib:

  A `BernoulliDistrib` object, from
  [`bernoulli_distrib()`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `y`. `mu` must lie in \\(0, 1)\\.

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

  A single positive integer, how many threads the kernel may use.
  Defaults to `1L`.

## Value

A named list of one numeric vector, `mu_mu`, of length
`max(length(y), length(mu))` and constant at \\-1/(\mu(1-\mu))\\.

## Notation

The **expected information** is
\\\mathbb{E}\[-\partial^2\ell/\partial\mu^2\]\\, the expectation of the
**observed information** under the model. The Bernoulli is a regular
family, so the second Bartlett identity holds and this equals the
variance of the score.

## See also

[`distrib_hessian.BernoulliDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.BernoulliDistrib.md)
for the observed quantity this is the expectation of,
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which inverts it, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- bernoulli_distrib()
th <- list(mu = 0.3)

# A single number, the reciprocal of the variance.
unique(distrib_expected_hessian(d, c(0, 1, 1), th)$mu_mu)
#> [1] -4.761905
-1 / (0.3 * 0.7)
#> [1] -4.761905

# Least information where the outcome is most uncertain.
vapply(c(0.5, 0.1, 0.01), function(p)
  -distrib_expected_hessian(d, 0, list(mu = p))$mu_mu, numeric(1))
#> [1]   4.00000  11.11111 101.01010

# On the canonical logit link the observed and the expected values agree
# exactly, at every observation.
rbind(observed = distrib_hessian(d, c(0, 1, 1), th, scale = "link")$mu_mu,
      expected = distrib_expected_hessian(d, c(0, 1, 1), th,
                                          scale = "link")$mu_mu)
#>           [,1]  [,2]  [,3]
#> observed -0.21 -0.21 -0.21
#> expected -0.21 -0.21 -0.21
```
