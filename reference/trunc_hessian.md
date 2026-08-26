# Observed Hessian of a Truncated Distribution

Computes \$\$\frac{\partial^2 \ell_T}{\partial\theta_i \partial\theta_j}
= H\_{ij}(y) - M\_{ij} + m_i m_j,\$\$ the subtracted part \\M\_{ij} -
m_i m_j\\ being the second derivative of \\\log Z\\. One of the shared
bodies, registered on both truncated classes.

## Usage

``` r
trunc_hessian(distrib, y, theta, scale = c("parameter", "link"), ...)
```

## Arguments

- distrib:

  A truncated distribution object, of either class.

- y:

  A numeric vector of observations.

- theta:

  A named list of the parent's parameters.

- scale:

  One of `"parameter"` (the default) or `"link"`, applied by the generic
  before dispatch.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of length `length(y)`, one component per
unordered pair of parameters, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## Details

Only \\H\_{ij}(y)\\ varies with the data; the correction is a function
of \\\theta\\ and the endpoints alone, so it is computed once and
recycled to the length of `y`.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`trunc_M()`](https://statmodels7.github.io/distributions7/reference/trunc_M.md)
and
[`trunc_score_mean()`](https://statmodels7.github.io/distributions7/reference/trunc_score_mean.md)
for the two corrections,
[`trunc_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/trunc_expected_hessian.md)
for the expectation, and
[`distrib_hessian.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.TruncatedContinuousDistrib.md)
and
[`distrib_hessian.TruncatedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.TruncatedDiscreteDistrib.md),
the two registrations.

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

y <- c(-0.5, 0, 1)
H <- distrib_hessian(tn, y, theta)
vapply(H, sum, numeric(1))
#>       mu_mu sigma_sigma    mu_sigma 
#>  -0.8721839   0.4749147   0.6374049 

# Against a numerical Hessian of the truncated log-likelihood.
ll <- function(v) {
  t2 <- as.list(v); names(t2) <- tn@params
  sum(distrib_pdf(tn, y, t2, log = TRUE))
}
Hn <- numDeriv::hessian(ll, unlist(theta))
ref <- vapply(distributions7:::hess_pairs(tn@params),
              function(q) Hn[match(q[1], tn@params), match(q[2], tn@params)],
              numeric(1))
max(abs(vapply(H, sum, numeric(1)) - ref))
#> [1] 6.119927e-11

# The correction is one number per component, not one per observation.
Hp <- distrib_hessian(gaussian1_distrib(), y, theta)
round(unlist(H) - unlist(Hp), 8)
#>       mu_mu1       mu_mu2       mu_mu3 sigma_sigma1 sigma_sigma2 sigma_sigma3 
#>   0.40371646   0.40371646   0.40371646   0.05220921   0.05220921   0.05220921 
#>    mu_sigma1    mu_sigma2    mu_sigma3 
#>   0.05814730   0.05814730   0.05814730 
```
