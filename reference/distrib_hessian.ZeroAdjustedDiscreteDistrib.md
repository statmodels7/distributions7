# Zero-Adjusted Discrete Observed Hessian

Computes the second derivatives of the hurdle log-mass in closed form.
The MIXED BLOCKS ARE IDENTICALLY ZERO, at every observation and every
parameter value, because the likelihood factorizes into a binary part in
\\\pi\\ and a positive-count part in \\\theta\\. The parent block is
zero at a zero and, at a positive observation, the parent's own Hessian
plus the truncation correction \$\$H\_{\mathrm{corr}} = \frac{\\1 -
f(0)\\\\f''(0) + f'(0)^2}{\\1 - f(0)\\^2},\$\$ the second derivative of
\\-\log\\1 - f(0)\\\\. The hurdle block is \\-1/\pi^2\\ at a zero and
\\-1/(1-\pi)^2\\ elsewhere.

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
  before dispatch.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors, one per unordered pair of parameters,
keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md).
Every mixed key holds a vector of exact zeros.

## Notation

\\f\\ is the parent's mass function, \\\pi\\ the probability of a zero,
\\s\\ the parent's score, \\H\\ its observed Hessian and \\\ell\\ the
log-mass of one observation.

## See also

[`distrib_gradient.ZeroAdjustedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ZeroAdjustedDiscreteDistrib.md)
for the first order,
[`distrib_expected_hessian.ZeroAdjustedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.ZeroAdjustedDiscreteDistrib.md)
for the expectation,
[`distrib_hessian.ZeroInflatedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ZeroInflatedDistrib.md),
whose mixed block is NOT zero, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- zero_adjusted(poisson_distrib())
theta <- list(mu = 3, za = 0.4)
set.seed(2)
y <- distrib_rng(d, 300, theta)

H <- distrib_hessian(d, y, theta)
vapply(H, sum, numeric(1))
#>       mu_mu       za_za       mu_za 
#>   -48.17993 -1298.61111     0.00000 

# Against a numerical Hessian of the log-likelihood.
ll <- function(v) {
  t2 <- as.list(v); names(t2) <- d@params
  sum(distrib_pdf(d, y, t2, log = TRUE))
}
Hn <- numDeriv::hessian(ll, unlist(theta))
ref <- vapply(distributions7:::hess_pairs(d@params),
              function(q) Hn[match(q[1], d@params), match(q[2], d@params)],
              numeric(1))
max(abs(vapply(H, sum, numeric(1)) - ref))
#> [1] 3.996092e-09

# The mixed block is exactly zero at every observation, where the
# zero-inflated wrapper's is not.
all(H$mu_za == 0)
#> [1] TRUE
zi <- zero_inflated(poisson_distrib())
unique(distrib_hessian(zi, 0:2, list(mu = 3, zi = 0.4))$mu_zi != 0)
#> [1]  TRUE FALSE
```
