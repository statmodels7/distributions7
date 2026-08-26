# Geometric Expected Hessian

Returns the expectation of the observed second derivative under the
model, in closed form and with no quadrature or simulation:
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{1}{\mu(1+\mu)},\$\$ which follows from \\\mathbb{E}\[Y\] = \mu\\
killing the data term. The Fisher information is \\1/(\mu(1+\mu))\\, the
reciprocal of the variance, so it is smaller than a Poisson's \\1/\mu\\
at the same mean: the extra dispersion is paid for in precision.

Because the value does not depend on the data, `approx` and `nsim` are
ignored. `y` is read only for its length.

## Arguments

- distrib:

  A `GeometricDistrib` object, from
  [`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md).

- y:

  A numeric vector of counts. Only its length is used.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `y`. `mu` must be strictly positive.

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
`max(length(y), length(mu))` and constant at \\-1/(\mu(1+\mu))\\.

## Notation

The **expected information** is
\\\mathbb{E}\[-\partial^2\ell/\partial\mu^2\]\\, the expectation of the
**observed information** under the model. The geometric is a regular
family, so the second Bartlett identity holds and this equals the
variance of the score.

## See also

[`distrib_hessian.GeometricDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GeometricDistrib.md)
for the observed quantity this is the expectation of,
[`distrib_expected_hessian.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.PoissonDistrib.md)
for the equidispersed comparison,
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which inverts it, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- geometric_distrib()
th <- list(mu = 3)

# A single number, the reciprocal of the variance.
unique(distrib_expected_hessian(d, c(0, 2, 7), th)$mu_mu)
#> [1] -0.08333333
-1 / (3 * (1 + 3))
#> [1] -0.08333333

# Less information than a Poisson of the same mean, the variance being
# larger.
c(geometric = -distrib_expected_hessian(d, 0, th)$mu_mu,
  poisson = -distrib_expected_hessian(poisson_distrib(), 0, th)$mu_mu)
#>  geometric    poisson 
#> 0.08333333 0.33333333 

# The observed value averages onto it over a large sample.
set.seed(5)
z <- distrib_rng(d, 2e5, th)
c(observed = mean(distrib_hessian(d, z, th)$mu_mu),
  expected = distrib_expected_hessian(d, 0, th)$mu_mu)
#>    observed    expected 
#> -0.08388410 -0.08333333 
```
