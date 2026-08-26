# Laplace Expected Hessian, Rate Parametrization

Returns the negative of the Fisher information, in closed form and with
no quadrature or simulation: \$\$-I(\mu) = -\lambda^2, \qquad
-I(\lambda) = -\dfrac{1}{\lambda^2}, \qquad -I(\mu, \lambda) = 0.\$\$

For \\\mu\\ this is the **variance of the score**. It is not the
expectation of
[`distrib_hessian.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Laplace2Distrib.md),
which is identically zero: the family has a kink at \\y = \mu\\ and the
second Bartlett identity fails there, exactly as in
[`distrib_expected_hessian.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.LaplaceDistrib.md).
With \\\partial\ell/\partial\mu = \lambda\\\mathrm{sign}(r)\\ the
variance is \\\lambda^2\\. The mixed entry vanishes because
\\\mathbb{E}\[\mathrm{sign}(r)\] = 0\\ by symmetry, and the rate entry
does equal its observed counterpart, that one carrying no data.

Because the values do not depend on the data, `approx` and `nsim` are
ignored. `y` is read only for its length.

## Arguments

- distrib:

  A `Laplace2Distrib` object, from
  [`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with components `mu` and `lambda`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `lambda` must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored, a closed form being available. Accepted so that the signature
  matches the generic's. As for
  [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md),
  a strategy averaging the observed second derivative would answer 0 in
  the location.

- nsim:

  Ignored, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of three numeric vectors, `mu_mu`, `lambda_lambda` and
`mu_lambda`, each of length `length(y)`.

## Notation

The **expected information** is the variance of the score. For a regular
family it also equals the expectation of the **observed information**,
by the second Bartlett identity; this family is not regular in \\\mu\\
and the two differ there.

## See also

[`distrib_hessian.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Laplace2Distrib.md),
whose `mu_mu` is zero;
[`distrib_expected_hessian.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.LaplaceDistrib.md)
for the same argument in the scale parametrization;
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which inverts this matrix; and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- laplace2_distrib()
th <- list(mu = 0.4, lambda = 2)

# -lambda^2 in the location, -1/lambda^2 in the rate, 0 mixed.
lapply(distrib_expected_hessian(d, c(-1.2, 0.3, 2.5), th), unique)
#> $mu_mu
#> [1] -4
#> 
#> $lambda_lambda
#> [1] -0.25
#> 
#> $mu_lambda
#> [1] 0
#> 
c(-2^2, -1 / 2^2)
#> [1] -4.00 -0.25

# The location entry is the variance of the score.
set.seed(12)
z <- distrib_rng(d, 1e5, th)
c(var_of_score = mean(distrib_gradient(d, z, th)$mu^2),
  information = -distrib_expected_hessian(d, 0, th)$mu_mu)
#> var_of_score  information 
#>            4            4 

# It agrees with the scale parametrization once the chain rule is applied:
# I(mu) is the same number, the location being shared.
-distrib_expected_hessian(laplace_distrib(), 0,
                          list(mu = 0.4, sigma = 1 / 2))$mu_mu
#> [1] 4
```
