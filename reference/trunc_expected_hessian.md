# Expected Hessian of a Truncated Distribution

Computes \$\$\mathbb{E}\left\[\frac{\partial^2 \ell_T}{\partial\theta_i
\partial\theta_j}\right\] = -\mathrm{Cov}\_T(s_i, s_j),\$\$ the second
Bartlett identity under the truncated law. The truncated-law mean of the
parent's Hessian, \\\mathbb{E}\_T\[H\_{ij}\]\\, cancels against the same
quantity in \\M\_{ij}\\, which enters the second derivative of \\\log
Z\\, and \\-\mathrm{Cov}\_T(s_i, s_j)\\ remains. The body is shared and
is registered on both truncated classes.

## Usage

``` r
trunc_expected_hessian(
  distrib,
  y,
  theta,
  scale = c("parameter", "link"),
  approx = c("opg", "bartlett", "integrate", "mc"),
  nsim = 10000,
  ...
)
```

## Arguments

- distrib:

  A truncated distribution object, of either class.

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list of the parent's parameters.

- scale:

  One of `"parameter"` (the default) or `"link"`. The generic applies
  the link scale to the value that the method returns.

- approx:

  Ignored. It is present so that the signature matches the generic's.

- nsim:

  Ignored, for the same reason.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of length `length(y)`, one component per
unordered pair of parameters, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md);
each vector is constant for scalar parameters.

## Details

For a parent without a scalar route this needs one quadrature per
component even where the parent has exact cdf derivatives. Those give
\\d\_{ij} Z / Z\\, which is \\\mathbb{E}\_T\[H\_{ij}\] +
\mathbb{E}\_T\[s_i s_j\]\\, and the covariance needs the second term on
its own. For a parent with a scalar route the sums come from
[`trunc_route_parts()`](https://statmodels7.github.io/distributions7/reference/trunc_route_parts.md).

The components do not depend on the data, so for a single parameter
vector every returned vector is constant; if a parameter varies by
observation, the components vary with it. The arguments `approx` and
`nsim` are accepted so that the signature matches the generic's, and the
method does not read them: the identity above is exact, and the
quantities on which it rests are computed by quadrature for every value
of `approx`.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, which both belong to the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`trunc_score_prod_mean()`](https://statmodels7.github.io/distributions7/reference/trunc_score_prod_mean.md)
and
[`trunc_score_mean()`](https://statmodels7.github.io/distributions7/reference/trunc_score_mean.md)
for the two pieces of the covariance,
[`trunc_hessian()`](https://statmodels7.github.io/distributions7/reference/trunc_hessian.md)
for the observed matrix, and
[`distrib_expected_hessian.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.TruncatedContinuousDistrib.md)
and
[`distrib_expected_hessian.TruncatedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.TruncatedDiscreteDistrib.md),
the two registrations.

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

EH <- distrib_expected_hessian(tn, 1, theta)
unlist(EH)
#>       mu_mu sigma_sigma    mu_sigma 
#> -0.29072798 -0.14498150 -0.07605547 

# The result equals minus the covariance of the parent's score under the
# truncated law.
m <- distributions7:::trunc_score_mean(tn, theta)
ES <- distributions7:::trunc_score_prod_mean(tn, theta)
c(mu_mu = -(ES$mu_mu - m$mu^2),
  sigma_sigma = -(ES$sigma_sigma - m$sigma^2),
  mu_sigma = -(ES$mu_sigma - m$mu * m$sigma))
#>       mu_mu sigma_sigma    mu_sigma 
#> -0.29072798 -0.14498150 -0.07605547 

# And what averaging the observed Hessian over the truncated law gives.
set.seed(1)
ys <- distrib_rng(tn, 200000, theta)
round(vapply(distrib_hessian(tn, ys, theta), mean, numeric(1)), 3)
#>       mu_mu sigma_sigma    mu_sigma 
#>      -0.291      -0.146      -0.075 
round(vapply(EH, function(z) z[1], numeric(1)), 3)
#>       mu_mu sigma_sigma    mu_sigma 
#>      -0.291      -0.145      -0.076 
```
