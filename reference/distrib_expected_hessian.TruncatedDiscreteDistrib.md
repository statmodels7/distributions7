# Truncated Analytical Expected Hessian (Discrete)

Computes \$\$\mathbb{E}\left\[\frac{\partial^2 \ell_T}
{\partial\theta_i\partial\theta_j}\right\] = -\mathrm{Cov}\_T(s_i,
s_j),\$\$ the covariance of the PARENT's score over the retained support
points. It is the second Bartlett identity for the truncated model, and
for a discrete parent the expectations behind it are exact sums over
that support.

## Arguments

- distrib:

  A `TruncatedDiscreteDistrib` object, from
  [`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list of the parent's parameters.

- scale:

  One of `"parameter"` (the default) or `"link"`, applied by the generic
  before dispatch.

- approx:

  Ignored. Present so that the signature matches the generic's.

- nsim:

  Ignored, for the same reason.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of length `length(y)`, one component per
unordered pair of parameters, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
each vector constant.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md),
[`trunc_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/trunc_expected_hessian.md)
for the shared body,
[`trunc_score_prod_mean()`](https://statmodels7.github.io/distributions7/reference/trunc_score_prod_mean.md)
for the moment it rests on, and
[`distrib_hessian.TruncatedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.TruncatedDiscreteDistrib.md)
for the observed matrix.

## Examples

``` r
ztp <- truncated(poisson_distrib(), lower = 1)
theta <- list(mu = 2)

EH <- distrib_expected_hessian(ztp, 1, theta)
EH$mu_mu
#> [1] -0.3972434

# Minus the variance of the parent's score over the retained support.
m <- distributions7:::trunc_score_mean(ztp, theta)$mu
ES <- distributions7:::trunc_score_prod_mean(ztp, theta)$mu_mu
-(ES - m^2)
#> [1] -0.3972434

# And what averaging the observed Hessian gives.
set.seed(9)
ys <- distrib_rng(ztp, 200000, theta)
c(sampled = mean(distrib_hessian(ztp, ys, theta)$mu_mu),
  closed = EH$mu_mu)
#>    sampled     closed 
#> -0.3967571 -0.3972434 
```
