# Truncated Score Mean by Quadrature

Computes \\m_i = \mathbb{E}\_T\[s_i\]\\ by integrating the parent's
score against the truncated law, one integral per parameter. It is the
route
[`trunc_score_mean()`](https://statmodels7.github.io/distributions7/reference/trunc_score_mean.md)
falls back on for a parent whose cdf derivatives are not exact, and is
the only route for a gamma or a beta, whose distribution function has no
elementary derivative in its shape.

## Usage

``` r
trunc_score_mean_quad(distrib, theta)
```

## Arguments

- distrib:

  A truncated distribution object, of either class.

- theta:

  A named list of the parent's parameters.

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

[`trunc_score_mean()`](https://statmodels7.github.io/distributions7/reference/trunc_score_mean.md),
which prefers the cdf route where it is exact, and
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md),
which performs the integration.

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

# The two routes agree; only the cost differs.
unlist(distributions7:::trunc_score_mean_quad(tn, theta))
#>          mu       sigma 
#>  0.08052166 -0.47667927 
unlist(distributions7:::trunc_score_mean(tn, theta))
#>          mu       sigma 
#>  0.08052166 -0.47667927 

# For a gamma parent this is the only route available.
tg <- truncated(gamma2_distrib(), lower = 0.5, upper = 5)
unlist(distributions7:::trunc_score_mean_quad(tg, list(mu = 2, sigma2 = 1)))
#>          mu      sigma2 
#>  0.07911249 -0.08650595 
```
