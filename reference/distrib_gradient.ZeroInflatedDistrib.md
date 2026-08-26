# Zero-Inflated Score

Computes the first derivatives of the zero-inflated log-mass in closed
form. For the parent's parameters the score is the parent's own,
weighted by the posterior probability that an observed zero came from
the parent, \$\$\frac{\partial \ell}{\partial \theta_i} = w\\ s_i(y),
\qquad w = \begin{cases} \dfrac{(1-\zeta)f(0)}{L_0} & y = 0 \\ 1 & y \>
0, \end{cases}\$\$ and for the inflation parameter \$\$\frac{\partial
\ell}{\partial \zeta} = \mathbb{I}(y = 0)\frac{1 - f(0)}{L_0} -
\mathbb{I}(y \> 0)\frac{1}{1 - \zeta}.\$\$

## Arguments

- distrib:

  A `ZeroInflatedDistrib` object, from
  [`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with the parent's parameters followed by `zi`.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch. `zi` rides a logit by default, so the two differ.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with one numeric vector per parameter, in `distrib@params`
order.

## Details

A positive observation carries no information about \\\zeta\\ beyond the
\\(1-\zeta)\\ factor, which is why its \\\zeta\\ component is the same
number at every such point. A zero carries all of it, and its weight
\\w\\ is what shares the credit between the two mechanisms.

## Notation

\\f\\ is the parent's mass function, \\\zeta\\ the probability of a
structural zero, \\L_0 = \zeta + (1-\zeta)f(0)\\ the inflated mass at
zero, \\w = (1-\zeta)f(0)/L_0\\ the posterior probability that an
observed zero came from the parent, \\s\\ the parent's score, \\H\\ its
observed Hessian and \\\ell\\ the log-mass of one observation.

## See also

[`distrib_hessian.ZeroInflatedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ZeroInflatedDistrib.md)
for the second order,
[`distrib_expected_hessian.ZeroInflatedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.ZeroInflatedDistrib.md)
for the information, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- zero_inflated(poisson_distrib())
theta <- list(mu = 3, zi = 0.25)
set.seed(2)
y <- distrib_rng(d, 200, theta)

g <- distrib_gradient(d, y, theta)
vapply(g, sum, numeric(1))
#>        mu        zi 
#> -3.299986 53.511184 

# Against a numerical derivative of the log-likelihood.
ll <- function(v) {
  t2 <- as.list(v); names(t2) <- d@params
  sum(distrib_pdf(d, y, t2, log = TRUE))
}
max(abs(vapply(g, sum, numeric(1)) - numDeriv::grad(ll, unlist(theta))))
#> [1] 8.314181e-09

# Every positive observation gives the same zi component, -1 / (1 - zi).
c(unique(round(g$zi[y > 0], 12)), theory = -1 / 0.75)
#>              theory 
#> -1.333333 -1.333333 

# And the parent's component is zero at a zero only in the limit: it is
# the parent's score there, weighted down by w.
c(weighted = g$mu[y == 0][1], unweighted = distrib_gradient(
    poisson_distrib(), 0, list(mu = 3))$mu)
#>   weighted unweighted 
#> -0.1299515 -1.0000000 
```
