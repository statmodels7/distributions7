# Zero-Adjusted Discrete Score

Computes the first derivatives of the hurdle log-mass in closed form.
The likelihood SEPARATES, so each parameter is informed by one half of
the data only. For the parent's parameters the score is zero at a zero
and the truncated score at a positive observation, \$\$\frac{\partial
\ell}{\partial \theta_i} = \mathbb{I}(y \> 0)\left\\s_i(y) +
\frac{f(0)}{1 - f(0)}\\s_i(0)\right\\,\$\$ the second term being the
derivative of the normalizing constant \\-\log\\1 - f(0)\\\\. For the
hurdle parameter, \$\$\frac{\partial \ell}{\partial \pi} =
\frac{\mathbb{I}(y = 0)}{\pi} - \frac{\mathbb{I}(y \> 0)}{1 - \pi},\$\$
which is the Bernoulli score of the indicator that the observation is
zero.

## Arguments

- distrib:

  A `ZeroAdjustedDiscreteDistrib` object, from
  [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with the parent's parameters followed by `za`.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch. `za` rides a logit by default.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with one numeric vector per parameter, in `distrib@params`
order.

## Notation

\\f\\ is the parent's mass function, \\\pi\\ the probability of a zero,
\\s\\ the parent's score, \\H\\ its observed Hessian and \\\ell\\ the
log-mass of one observation.

## See also

[`distrib_hessian.ZeroAdjustedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ZeroAdjustedDiscreteDistrib.md)
for the second order,
[`distrib_gradient.ZeroInflatedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ZeroInflatedDistrib.md),
where the two halves do NOT separate, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- zero_adjusted(poisson_distrib())
theta <- list(mu = 3, za = 0.4)
set.seed(2)
y <- distrib_rng(d, 300, theta)

g <- distrib_gradient(d, y, theta)
vapply(g, sum, numeric(1))
#>        mu        za 
#> -2.697686 58.333333 

# Against a numerical derivative of the log-likelihood.
ll <- function(v) {
  t2 <- as.list(v); names(t2) <- d@params
  sum(distrib_pdf(d, y, t2, log = TRUE))
}
max(abs(vapply(g, sum, numeric(1)) - numDeriv::grad(ll, unlist(theta))))
#> [1] 2.340322e-09

# A zero tells the parent's parameters nothing at all, the two halves of
# the likelihood being separate.
unique(g$mu[y == 0])
#> [1] 0

# And the za component is the Bernoulli score of the zero indicator, so it
# takes exactly two values.
unique(round(g$za, 10))
#> [1]  2.500000 -1.666667
c(1 / 0.4, -1 / 0.6)
#> [1]  2.500000 -1.666667
```
