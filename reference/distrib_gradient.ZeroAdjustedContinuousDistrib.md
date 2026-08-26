# Zero-Adjusted Continuous Score

Computes the first derivatives of the mixed log-density in closed form.
The likelihood separates COMPLETELY here, with no truncation term to
correct for: the parent's parameters take the parent's own score at a
non-zero observation and zero at the atom, and \$\$\frac{\partial
\ell}{\partial \pi} = \frac{\mathbb{I}(y = 0)}{\pi} - \frac{\mathbb{I}(y
\ne 0)}{1 - \pi},\$\$ the Bernoulli score of the indicator that the
observation is zero.

## Arguments

- distrib:

  A `ZeroAdjustedContinuousDistrib` object, from
  [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with the parent's parameters followed by `za`.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with one numeric vector per parameter, in `distrib@params`
order.

## Details

Compare the discrete branch, whose parent block carries the truncation
correction \\f(0)s(0)/\\1 - f(0)\\\\. There is none here because a
continuous parent places no mass at zero, so nothing is removed from it.

## Notation

\\f_W\\ is the parent's density, \\\pi\\ the probability of the atom at
zero, \\f_Y\\ the mixed density and \\\ell\\ the log-density of one
observation.

## See also

[`distrib_hessian.ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ZeroAdjustedContinuousDistrib.md)
for the second order,
[`distrib_gradient.ZeroAdjustedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ZeroAdjustedDiscreteDistrib.md),
which carries the truncation term, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- zero_adjusted(gaussian1_distrib())
theta <- list(mu = 1, sigma = 2, za = 0.3)
set.seed(4)
y <- distrib_rng(d, 300, theta)

g <- distrib_gradient(d, y, theta)
vapply(g, sum, numeric(1))
#>        mu     sigma        za 
#> -6.757865 -4.787530 42.857143 

# Against a numerical derivative of the log-likelihood.
ll <- function(v) {
  t2 <- as.list(v); names(t2) <- d@params
  sum(distrib_pdf(d, y, t2, log = TRUE))
}
max(abs(vapply(g, sum, numeric(1)) - numDeriv::grad(ll, unlist(theta))))
#> [1] 4.47217e-09

# Away from the atom it is the parent's own score, with no correction.
nz <- which(y != 0)[1]
c(mixed = g$mu[nz],
  parent = distrib_gradient(gaussian1_distrib(), y[nz],
                            theta[c("mu", "sigma")])$mu)
#>     mixed    parent 
#> 0.4384923 0.4384923 

# And an observation at the atom says nothing about the parent.
unique(g$mu[y == 0])
#> [1] 0
```
