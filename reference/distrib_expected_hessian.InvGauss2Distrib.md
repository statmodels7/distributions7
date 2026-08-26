# Inverse Gaussian Expected Hessian in Mean and Shape

Returns the expectation of the observed Hessian under the model, in
closed form and with no quadrature or simulation. Every second
derivative is at most linear in the response, so the expectations need
only \\\mathbb{E}\[Y\] = \mu\\:
\$\$\mathbb{E}\left\[\ell^{(\mu\mu)}\right\] = -\dfrac{\lambda}{\mu^3},
\qquad \mathbb{E}\left\[\ell^{(\mu\lambda)}\right\] = 0, \qquad
\mathbb{E}\left\[\ell^{(\lambda\lambda)}\right\] =
-\dfrac{1}{2\lambda^2}.\$\$ The pure shape entry is the observed value
itself, being free of the data.

Both diagonal entries are negative at every parameter setting, so the
information is positive definite everywhere. The zero off-diagonal says
the mean and the shape are orthogonal, so the mean equation can be
fitted with the shape held at any value without biasing it.

Because the values do not depend on the data, `approx` and `nsim` are
ignored. `y` is read only for its length.

## Arguments

- distrib:

  An `InvGauss2Distrib` object, from
  [`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with components `mu` and `lambda`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

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

A named list of three numeric vectors, `mu_mu`, `mu_lambda` and
`lambda_lambda`, in that order, each of length
`max(length(y), length(mu), length(lambda))` and constant within itself
when the parameters are.

## Notation

The **expected information** is
\\\mathbb{E}\[-\partial^2\ell/\partial\theta\\\partial\theta^\top\]\\,
the expectation of the **observed information** under the model. The
inverse Gaussian is a regular family, so the second Bartlett identity
holds and this equals the variance of the score. \\\lambda\\ names this
family's shape parameter.

## See also

[`distrib_hessian.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.InvGauss2Distrib.md)
for the observed quantity this is the expectation of,
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which inverts it at each step, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- invgauss2_distrib()
th <- list(mu = 2, lambda = 3)

# The three constants, one value per observation.
lapply(distrib_expected_hessian(d, c(1, 2, 3), th), unique)
#> $mu_mu
#> [1] -0.375
#> 
#> $mu_lambda
#> [1] 0
#> 
#> $lambda_lambda
#> [1] -0.05555556
#> 
c(-3 / 2^3, 0, -1 / (2 * 3^2))
#> [1] -0.37500000  0.00000000 -0.05555556

# The pure shape entry equals the observed one at every observation.
identical(distrib_expected_hessian(d, c(1, 2, 3), th)$lambda_lambda,
          distrib_hessian(d, c(1, 2, 3), th)$lambda_lambda)
#> [1] TRUE

# The observed Hessian averages onto them over a large sample.
set.seed(11)
z <- distrib_rng(d, 2e4, th)
rbind(observed = vapply(distrib_hessian(d, z, th), mean, numeric(1)),
      expected = vapply(distrib_expected_hessian(d, z, th),
                        function(v) v[1], numeric(1)))
#>               mu_mu   mu_lambda lambda_lambda
#> observed -0.3827783 0.001728511   -0.05555556
#> expected -0.3750000 0.000000000   -0.05555556
```
