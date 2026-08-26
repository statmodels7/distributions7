# Truncated Analytical Observed Hessian (Discrete)

Computes \$\$\frac{\partial^2 \ell_T}{\partial\theta_i\partial\theta_j}
= H\_{ij}(y) - M\_{ij} + m_i m_j, \qquad M\_{ij} =
\mathbb{E}\_T\[H\_{ij} + s_i s_j\],\$\$ the subtracted part being the
second derivative of \\\log Z\\. For a discrete parent \\M\_{ij}\\ comes
from the parent's cdf derivatives, which are exact finite sums, with the
mass at the retained lower endpoint added back.

## Arguments

- distrib:

  A `TruncatedDiscreteDistrib` object, from
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
[`distrib_expected_hessian.TruncatedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.TruncatedDiscreteDistrib.md)
for the expectation.

## Examples

``` r
ztp <- truncated(poisson_distrib(), lower = 1)
theta <- list(mu = 2)

y <- 1:4
H <- distrib_hessian(ztp, y, theta)
H$mu_mu
#> [1] -0.06898458 -0.31898458 -0.56898458 -0.81898458

# Against a numerical Hessian of the truncated log-likelihood.
ll <- function(v) {
  t2 <- as.list(v); names(t2) <- ztp@params
  sum(distrib_pdf(ztp, y, t2, log = TRUE))
}
abs(sum(H$mu_mu) - numDeriv::hessian(ll, unlist(theta))[1, 1])
#> [1] 1.200817e-12
```
