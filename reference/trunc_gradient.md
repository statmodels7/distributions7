# Score of a Truncated Distribution

Computes the parent's score recentered by its truncated mean,
\$\$\frac{\partial \ell_T}{\partial \theta_i} = s_i(y) - m_i, \qquad m_i
= \mathbb{E}\_T\[s_i\],\$\$ the subtraction being the derivative of
\\-\log Z\\. One of the shared bodies, registered on both truncated
classes.

## Usage

``` r
trunc_gradient(distrib, y, theta, scale = c("parameter", "link"), ...)
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

A named list of numeric vectors, one component per parameter, in
`distrib@params` order.

## Details

Truncation adds no parameter, the endpoints being constants, so the
returned list has exactly the parent's components. What it adds is the
\\\theta\\-dependent normalizing constant, and the recentering is that
constant's whole contribution at first order. The support does not
depend on \\\theta\\, which licenses differentiating \\Z\\ under the
integral sign and keeps truncation at fixed points a regular problem.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`trunc_score_mean()`](https://statmodels7.github.io/distributions7/reference/trunc_score_mean.md)
for \\m_i\\,
[`trunc_hessian()`](https://statmodels7.github.io/distributions7/reference/trunc_hessian.md)
for the second order, and
[`distrib_gradient.TruncatedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.TruncatedContinuousDistrib.md)
and
[`distrib_gradient.TruncatedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.TruncatedDiscreteDistrib.md),
the two registrations.

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

# It is the parent's score, shifted by one number per parameter.
gp <- distrib_gradient(gaussian1_distrib(), y, theta)
round(unlist(g) - unlist(gp), 8)
#>         mu1         mu2         mu3      sigma1      sigma2      sigma3 
#> -0.08052166 -0.08052166 -0.08052166  0.47667927  0.47667927  0.47667927 
```
