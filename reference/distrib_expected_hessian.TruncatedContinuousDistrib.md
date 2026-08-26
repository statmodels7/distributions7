# Truncated Analytical Expected Hessian (Continuous)

Computes \$\$\mathbb{E}\left\[\frac{\partial^2 \ell_T}
{\partial\theta_i\partial\theta_j}\right\] = -\mathrm{Cov}\_T(s_i,
s_j),\$\$ the covariance of the PARENT's score under the truncated law.
This is the second Bartlett identity for the truncated model: the
parent's expected Hessian cancels exactly against the term it
contributes to the second derivative of \\\log Z\\.

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object, from
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

## Details

One quadrature per component is needed even where the parent has exact
cdf derivatives, since those give \\\mathbb{E}\_T\[H\_{ij}\] +
\mathbb{E}\_T\[s_i s_j\]\\ and the covariance needs the second term on
its own. No component depends on the data, so every returned vector is
constant.

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
[`distrib_hessian.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.TruncatedContinuousDistrib.md)
for the observed matrix.

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

EH <- distrib_expected_hessian(tn, 1, theta)
unlist(EH)
#>       mu_mu sigma_sigma    mu_sigma 
#> -0.29072798 -0.14498150 -0.07605547 

# Minus the covariance of the parent's score under the truncated law.
m <- distributions7:::trunc_score_mean(tn, theta)
ES <- distributions7:::trunc_score_prod_mean(tn, theta)
c(mu_mu = -(ES$mu_mu - m$mu^2),
  sigma_sigma = -(ES$sigma_sigma - m$sigma^2),
  mu_sigma = -(ES$mu_sigma - m$mu * m$sigma))
#>       mu_mu sigma_sigma    mu_sigma 
#> -0.29072798 -0.14498150 -0.07605547 

# And what averaging the observed Hessian gives.
set.seed(1)
ys <- distrib_rng(tn, 200000, theta)
round(vapply(distrib_hessian(tn, ys, theta), mean, numeric(1)), 3)
#>       mu_mu sigma_sigma    mu_sigma 
#>      -0.291      -0.146      -0.075 
```
