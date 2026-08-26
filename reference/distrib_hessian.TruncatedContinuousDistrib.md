# Truncated Analytical Observed Hessian (Continuous)

Computes \$\$\frac{\partial^2 \ell_T}{\partial\theta_i\partial\theta_j}
= H\_{ij}(y) - M\_{ij} + m_i m_j, \qquad M\_{ij} =
\mathbb{E}\_T\[H\_{ij} + s_i s_j\],\$\$ the subtracted part being the
second derivative of \\\log Z\\. Only \\H\_{ij}(y)\\ varies with the
data; the correction is one number per component, recycled to the length
of `y`.

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object, from
  [`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md).

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

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md),
[`trunc_hessian()`](https://statmodels7.github.io/distributions7/reference/trunc_hessian.md)
for the shared body,
[`trunc_M()`](https://statmodels7.github.io/distributions7/reference/trunc_M.md)
for \\M\_{ij}\\, and
[`distrib_expected_hessian.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.TruncatedContinuousDistrib.md)
for the expectation.

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

# The correction is constant across observations.
round(unlist(H) - unlist(distrib_hessian(gaussian1_distrib(), y, theta)), 8)
#>       mu_mu1       mu_mu2       mu_mu3 sigma_sigma1 sigma_sigma2 sigma_sigma3 
#>   0.40371646   0.40371646   0.40371646   0.05220921   0.05220921   0.05220921 
#>    mu_sigma1    mu_sigma2    mu_sigma3 
#>   0.05814730   0.05814730   0.05814730 
```
