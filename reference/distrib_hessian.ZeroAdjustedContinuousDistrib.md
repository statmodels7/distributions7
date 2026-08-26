# Zero-Adjusted Continuous Observed Hessian

Computes the second derivatives of the mixed log-density in closed form.
The MIXED BLOCKS ARE EXACTLY ZERO, the likelihood factorizing into a
binary part and a continuous part. The parent block is the parent's own
observed Hessian at a non-zero observation and zero at the atom, with no
truncation correction; the hurdle block is \\-1/\pi^2\\ at the atom and
\\-1/(1-\pi)^2\\ elsewhere.

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

A named list of numeric vectors, one per unordered pair of parameters,
keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md).
Every key pairing a parent parameter with `za` holds a vector of exact
zeros.

## Notation

\\f_W\\ is the parent's density, \\\pi\\ the probability of the atom at
zero, \\f_Y\\ the mixed density and \\\ell\\ the log-density of one
observation.

## See also

[`distrib_gradient.ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.ZeroAdjustedContinuousDistrib.md)
for the first order,
[`distrib_expected_hessian.ZeroAdjustedContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.ZeroAdjustedContinuousDistrib.md)
for the expectation, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- zero_adjusted(gaussian1_distrib())
theta <- list(mu = 1, sigma = 2, za = 0.3)
set.seed(4)
y <- distrib_rng(d, 300, theta)

H <- distrib_hessian(d, y, theta)
vapply(H, sum, numeric(1))
#>        mu_mu  sigma_sigma        za_za     mu_sigma        mu_za     sigma_za 
#>   -50.250000   -93.318705 -1510.204082     6.757865     0.000000     0.000000 

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
#> [1] 4.180265e-09

# Every mixed block is exactly zero.
c(mu_za = all(H$mu_za == 0), sigma_za = all(H$sigma_za == 0))
#>    mu_za sigma_za 
#>     TRUE     TRUE 
```
