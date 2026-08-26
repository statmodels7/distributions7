# Truncated Analytical Gradient (Continuous)

Computes \$\$\frac{\partial \ell_T}{\partial\theta_i} = s_i(y) - m_i,
\qquad m_i = \mathbb{E}\_T\[s_i\],\$\$ the parent's score recentered at
its mean over the truncated support. The subtraction is the derivative
of \\-\log Z\\, and it gives the truncated score mean zero, which the
parent's own score does not have here.

## Arguments

- distrib:

  A `TruncatedContinuousDistrib` object, from
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
[`distrib_hessian.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.TruncatedContinuousDistrib.md)
for the second order.

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

y <- c(-0.5, 0, 1)
g <- distrib_gradient(tn, y, theta)
vapply(g, sum, numeric(1))
#>         mu      sigma 
#> -0.5193428 -0.3639437 

# Against a numerical derivative of the truncated log-likelihood.
ll <- function(v) {
  t2 <- as.list(v); names(t2) <- tn@params
  sum(distrib_pdf(tn, y, t2, log = TRUE))
}
max(abs(vapply(g, sum, numeric(1)) - numDeriv::grad(ll, unlist(theta))))
#> [1] 3.190437e-11

# The truncated score has mean zero; the parent's, over the same interval,
# does not.
set.seed(1)
ys <- distrib_rng(tn, 50000, theta)
round(vapply(distrib_gradient(tn, ys, theta), mean, numeric(1)), 3)
#>    mu sigma 
#> 0.001 0.000 
round(vapply(distrib_gradient(gaussian1_distrib(), ys, theta), mean,
             numeric(1)), 3)
#>     mu  sigma 
#>  0.081 -0.476 
```
