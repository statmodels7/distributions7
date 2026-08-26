# Mean of the Parent's Score Under the Truncated Law

Computes \\m_i = \mathbb{E}\_T\[s_i\]\\, the quantity that recenters the
parent's score: the truncated score is \\s_i(y) - m_i\\. It is \\d_i Z /
Z\\, because \\d_i Z = \int_T f s_i\\, so it comes from
[`trunc_mass_derivs()`](https://statmodels7.github.io/distributions7/reference/trunc_mass_derivs.md)
where those are exact and from
[`trunc_score_mean_quad()`](https://statmodels7.github.io/distributions7/reference/trunc_score_mean_quad.md)
otherwise.

## Usage

``` r
trunc_score_mean(distrib, theta)
```

## Arguments

- distrib:

  A truncated distribution object, of either class.

- theta:

  A named list of the parent's parameters.

## Value

A named list of numeric vectors, one component per parameter, in
`distrib@params` order.

## Details

\\m_i\\ is a function of \\\theta\\ and of the endpoints alone, not of
the data, so one evaluation serves a whole sample. It is also what makes
derivatives of a truncated distribution dearer than the parent's: the
quadrature route costs one integral per parameter.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`trunc_gradient()`](https://statmodels7.github.io/distributions7/reference/trunc_gradient.md),
which subtracts it,
[`trunc_M()`](https://statmodels7.github.io/distributions7/reference/trunc_M.md)
for the second-order counterpart, and
[`trunc_score_mean_quad()`](https://statmodels7.github.io/distributions7/reference/trunc_score_mean_quad.md)
for the fallback.

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

m <- distributions7:::trunc_score_mean(tn, theta)
unlist(m)
#>          mu       sigma 
#>  0.08052166 -0.47667927 

# The truncated score is the parent's, recentered by exactly this.
y <- c(-0.5, 0, 1)
g <- distrib_gradient(tn, y, theta)
gp <- distrib_gradient(gaussian1_distrib(), y, theta)
max(abs(unlist(g) - (unlist(gp) - rep(unlist(m), each = length(y)))))
#> [1] 0

# It is also what makes the truncated score have mean zero, which the
# untruncated score does not have over the interval.
set.seed(1)
ys <- distrib_rng(tn, 50000, theta)
round(vapply(distrib_gradient(tn, ys, theta), mean, numeric(1)), 3)
#>    mu sigma 
#> 0.001 0.000 
```
