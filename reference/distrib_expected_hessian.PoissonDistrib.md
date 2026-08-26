# Poisson Expected Hessian

Returns the expectation of the observed second derivative under the
model, in closed form and with no quadrature or simulation:
\$\$\mathbb{E}\left\[\dfrac{\partial^2 \ell}{\partial \mu^2}\right\] =
-\dfrac{1}{\mu},\$\$ which follows from \\\mathbb{E}\[Y\] = \mu\\ in
\\-y/\mu^2\\. The Fisher information for one observation is \\1/\mu\\,
so the asymptotic standard error of \\\hat\mu\\ is \\\sqrt{\mu/n}\\, the
standard error of a sample mean whose variance is \\\mu\\.

On the **link** scale with the default logarithm the value is \\-\mu\\,
and so is the observed Hessian: the log is the canonical link, the
observed curvature there carries no data, and the two coincide exactly.
Fisher scoring and Newton's method therefore take the same step on a
Poisson model with a log link, which is why the two are not
distinguished in the classical generalized-linear-model literature for
this case.

Because the value does not depend on the data, `approx` and `nsim` are
ignored. `y` is read only for its length.

## Arguments

- distrib:

  A `PoissonDistrib` object, from
  [`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md).

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
`max(length(y), length(mu))` and constant at \\-1/\mu\\.

## Notation

The **expected information** is
\\\mathbb{E}\[-\partial^2\ell/\partial\mu^2\]\\, the expectation of the
**observed information** under the model. The Poisson is a regular
family, so the second Bartlett identity holds and this equals the
variance of the score.

## See also

[`distrib_hessian.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.PoissonDistrib.md)
for the observed quantity this is the expectation of,
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which inverts it, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- poisson_distrib()
th <- list(mu = 3)

# A single number, -1/mu.
unique(distrib_expected_hessian(d, c(0, 2, 7), th)$mu_mu)
#> [1] -0.3333333
-1 / 3
#> [1] -0.3333333

# On the canonical log link the observed and the expected values agree
# exactly, at every observation.
rbind(observed = distrib_hessian(d, c(0, 2, 7), th, scale = "link")$mu_mu,
      expected = distrib_expected_hessian(d, c(0, 2, 7), th,
                                          scale = "link")$mu_mu)
#>          [,1] [,2] [,3]
#> observed   -3   -3   -3
#> expected   -3   -3   -3

# It is the variance of the score, this family being regular.
set.seed(4)
z <- distrib_rng(d, 2e5, th)
c(var_of_score = mean(distrib_gradient(d, z, th)$mu^2),
  information = -distrib_expected_hessian(d, 0, th)$mu_mu)
#> var_of_score  information 
#>    0.3326167    0.3333333 
```
