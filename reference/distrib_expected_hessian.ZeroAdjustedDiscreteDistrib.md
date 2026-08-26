# Zero-Adjusted Discrete Expected Information

Computes the expectation of the observed Hessian in closed form. The
hurdle block is the Bernoulli information \\-1/\\\pi(1-\pi)\\\\, the
mixed blocks are exactly zero, and the parent block is \\(1-\pi)\\ times
the truncated parent's own expected Hessian: a positive observation
arrives with probability \\1-\pi\\ and carries the whole of what
\\\theta\\ is estimated from. No component depends on the data, so every
returned vector is constant.

## Arguments

- distrib:

  A `ZeroAdjustedDiscreteDistrib` object, from
  [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with the parent's parameters followed by `za`.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch.

- approx:

  Ignored: the expectation is exact. Present so that the signature
  matches the generic's.

- nsim:

  Ignored, for the same reason.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of length `length(y)`, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
each vector constant.

## Details

The block diagonal is why a hurdle model is cheap to fit: the two halves
can be estimated independently, and the information matrix inverts
blockwise. `approx` and `nsim` are accepted so that the signature
matches the generic's and neither is read; the expectation here is
exact.

## Notation

\\f\\ is the parent's mass function, \\\pi\\ the probability of a zero,
\\s\\ the parent's score, \\H\\ its observed Hessian and \\\ell\\ the
log-mass of one observation.

## See also

[`distrib_hessian.ZeroAdjustedDiscreteDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ZeroAdjustedDiscreteDistrib.md)
for the observed matrix,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
whose Fisher scoring inverts this, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- zero_adjusted(poisson_distrib())
theta <- list(mu = 3, za = 0.4)
set.seed(2)
y <- distrib_rng(d, 300, theta)

EH <- distrib_expected_hessian(d, y, theta)
vapply(EH, function(z) z[1], numeric(1))
#>      mu_mu      za_za      mu_za 
#> -0.1773945 -4.1666667  0.0000000 

# The hurdle block is the Bernoulli information, and the mixed one is zero.
c(reported = EH$za_za[1], bernoulli = -1 / (0.4 * 0.6),
  mixed = EH$mu_za[1])
#>  reported bernoulli     mixed 
#> -4.166667 -4.166667  0.000000 

# It is what summing the observed Hessian against the mass function gives.
sup <- 0:400
m <- distrib_pdf(d, sup, theta)
Hs <- distrib_hessian(d, sup, theta)
rbind(summed = vapply(Hs, function(z) sum(z * m), numeric(1)),
      closed = vapply(EH, function(z) z[1], numeric(1)))
#>             mu_mu     za_za mu_za
#> summed -0.1773945 -4.166667     0
#> closed -0.1773945 -4.166667     0
```
