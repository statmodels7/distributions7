# Truncated Analytical Gradient (Discrete)

Computes \$\$\frac{\partial \ell_T}{\partial\theta_i} = s_i(y) - m_i,
\qquad m_i = \mathbb{E}\_T\[s_i\],\$\$ the parent's score recentered at
its mean over the retained support points. For a discrete parent \\m_i\\
is an exact finite sum, not a quadrature, the parent's cdf derivatives
being exact sums themselves.

## Arguments

- distrib:

  A `TruncatedDiscreteDistrib` object, from
  [`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md).

- y:

  A numeric vector of observations.

- theta:

  A named list of the parent's parameters. Truncation adds none, so the
  returned list has exactly the parent's components.

- scale:

  One of `"parameter"` (the default) or `"link"`, applied by the generic
  before dispatch.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors, one component per parameter, in
`distrib@params` order.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md),
[`trunc_gradient()`](https://statmodels7.github.io/distributions7/reference/trunc_gradient.md)
for the shared body,
[`trunc_score_mean()`](https://statmodels7.github.io/distributions7/reference/trunc_score_mean.md)
for \\m_i\\, and
[`distrib_hessian.TruncatedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.TruncatedDiscreteDistrib.md)
for the second order.

## Examples

``` r
ztp <- truncated(poisson_distrib(), lower = 1)
theta <- list(mu = 2)

y <- 1:4
g <- distrib_gradient(ztp, y, theta)
g$mu
#> [1] -0.6565176 -0.1565176  0.3434824  0.8434824

# The parent's score, shifted by one number.
distrib_gradient(poisson_distrib(), y, theta)$mu
#> [1] -0.5  0.0  0.5  1.0
distributions7:::trunc_score_mean(ztp, theta)$mu
#> [1] 0.1565176

# Against a numerical derivative of the truncated log-likelihood.
ll <- function(v) {
  t2 <- as.list(v); names(t2) <- ztp@params
  sum(distrib_pdf(ztp, y, t2, log = TRUE))
}
abs(sum(g$mu) - numDeriv::grad(ll, unlist(theta)))
#> [1] 2.370032e-11
```
