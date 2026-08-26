# Exponential Expected Hessian

Returns the expectation of the observed second derivative under the
model, in closed form and with no quadrature or simulation:
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{1}{\mu^2},\$\$ which follows from \\\mathbb{E}\[Y\] = \mu\\
substituted into \\(\mu - 2y)/\mu^3\\. The Fisher information for one
observation is \\1/\mu^2\\, so the asymptotic standard error of
\\\hat\mu\\ is \\\mu/\sqrt n\\, which is the standard error of a sample
mean of variance \\\mu^2\\.

Because the value does not depend on the data, `approx` and `nsim` are
ignored. `y` is read only for its length.

## Arguments

- distrib:

  An `ExponentialDistrib` object, from
  [`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

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
`max(length(y), length(mu))` and constant at \\-1/\mu^2\\.

## Notation

The **expected information** is
\\\mathbb{E}\[-\partial^2\ell/\partial\mu^2\]\\, the expectation of the
**observed information** under the model. The exponential is a regular
family, so the second Bartlett identity holds and this equals the
variance of the score.

## See also

[`distrib_hessian.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ExponentialDistrib.md)
for the observed quantity this is the expectation of,
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which inverts it, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- exponential_distrib()
th <- list(mu = 2)

# A single number, -1/mu^2.
unique(distrib_expected_hessian(d, c(0.3, 1.1, 4.0), th)$mu_mu)
#> [1] -0.25
-1 / 2^2
#> [1] -0.25

# The observed value averages onto it over a large sample.
set.seed(21)
z <- distrib_rng(d, 2e5, th)
c(observed = mean(distrib_hessian(d, z, th)$mu_mu),
  expected = distrib_expected_hessian(d, 0, th)$mu_mu)
#>   observed   expected 
#> -0.2498915 -0.2500000 

# It is the variance of the score, this family being regular.
c(var_of_score = mean(distrib_gradient(d, z, th)$mu^2),
  information = -distrib_expected_hessian(d, 0, th)$mu_mu)
#> var_of_score  information 
#>    0.2512417    0.2500000 
```
