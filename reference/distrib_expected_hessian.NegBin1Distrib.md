# NB1 Expected Hessian

Returns the expectation of the observed Hessian under the model. Every
term carrying \\P = \psi(y+r) - \psi(r) - \log(1+\theta)\\ drops out,
its expectation vanishing by the first Bartlett identity, and what
remains needs only \\\mathbb{E}\[\psi'(Y+r)\]\\. That has no closed
form: it is summed against the exact mass out to a far-tail quantile, so
this is a truncated exact sum and not a quadrature or a simulation.
`approx` and `nsim` are ignored, and `y` is read only for its length.

**The mixed entry does not vanish**, so the mean and the dispersion are
not orthogonal in this family. Measured at four settings it is 0.0172,
0.0364, 0.0108 and 0.0157; in
[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
the same entry is exactly zero. That is a difference between the two
negative binomials rather than a difference of parametrization.

## Arguments

- distrib:

  A `NegBin1Distrib` object, from
  [`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md).

- y:

  A numeric vector of counts. Only its length is used.

- theta:

  A named list with components `mu` and `theta`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored. The surviving expectation is a truncated exact sum. Accepted
  so that the signature matches the generic's.

- nsim:

  Ignored, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count says. Defaults to `1L`.

## Value

A named list of three numeric vectors, `mu_mu`, `mu_theta` and
`theta_theta`, in that order, each of length
`max(length(y), length(mu), length(theta))` and constant within itself
when the parameters are.

## A caveat at small theta

The dispersion entry inherits the cancellation of
[`distrib_gradient.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.NegBin1Distrib.md)
and **is not rewritten to remove it**. The chain rule through \\r =
\mu/\theta\\ divides by \\\theta^2\\ and \\\theta^4\\, so what is left
of the digits runs out early. Measured at \\\mu = 4\\: the entry reads
\\-0.489\\ at \\\theta = 10^{-2}\\ and \\-0.500\\ at \\10^{-4}\\, then
\\+2.1\times 10^{2}\\ at \\10^{-6}\\ and \\+2.9\times 10^{8}\\ at
\\10^{-8}\\. The sign is impossible for an expected second derivative,
and the matrix is indefinite there, its determinant turning negative. A
Fisher scoring step taken in that regime is not reliable, and a nearly
equidispersed sample drives a fit into it.

## Notation

The **expected information** is
\\\mathbb{E}\[-\partial^2\ell/\partial\theta\\\partial\theta^\top\]\\,
the expectation of the **observed information** under the model. The
first Bartlett identity is \\\mathbb{E}\[\partial\ell/\partial\theta\] =
0\\.

## See also

[`distrib_hessian.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.NegBin1Distrib.md)
for the observed quantity this is the expectation of,
[`distrib_expected_hessian.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.NegBin2Distrib.md)
for the quadratic-variance family, whose mixed entry is exactly zero,
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md),
which inverts it at each step, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- negbin1_distrib()
th <- list(mu = 4, theta = 4)
e <- distrib_expected_hessian(d, c(0, 2, 6), th)
lapply(e, unique)
#> $mu_mu
#> [1] -0.06717466
#> 
#> $mu_theta
#> [1] 0.01717466
#> 
#> $theta_theta
#> [1] -0.01717466
#> 

# Negative definite at a moderate dispersion.
M <- matrix(c(e$mu_mu[1], e$mu_theta[1], e$mu_theta[1], e$theta_theta[1]), 2)
eigen(M, only.values = TRUE)$values
#> [1] -0.01184367 -0.07250565

# The mixed entry is not zero, so the mean and the dispersion are not
# orthogonal here; in the quadratic-variance family it is exactly zero.
c(nb1 = e$mu_theta[1],
  nb2 = distrib_expected_hessian(negbin2_distrib(), 0, th)$mu_theta)
#>        nb1        nb2 
#> 0.01717466 0.00000000 

# The dispersion entry loses its digits as theta goes to zero, and turns
# positive, which an expected second derivative cannot be.
vapply(c(1e-2, 1e-4, 1e-6, 1e-8),
       function(t) distrib_expected_hessian(d, 0,
                     list(mu = 4, theta = t))$theta_theta,
       numeric(1))
#> [1] -4.889418e-01 -4.999694e-01  2.113984e+02  2.890546e+08
```
