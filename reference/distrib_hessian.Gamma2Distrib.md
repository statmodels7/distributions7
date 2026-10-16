# Gamma Observed Hessian in Mean and Variance

Computes the three distinct second derivatives of the gamma log-density
with respect to \\\mu\\ and \\\sigma^2\\, one value per observation, in
closed form, by differentiating the score of
[`distrib_gradient.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gamma2Distrib.md)
once more. Every component carries the trigamma function of \\\alpha =
\mu^2/\sigma^2\\, because both parameters move the shape.

The curvature in \\\mu\\ is not negative at every data point, and the
same holds for the curvature in \\\sigma^2\\: the observed information
of a gamma is not positive definite everywhere. Its expectation is
positive definite; see
[`distrib_expected_hessian.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gamma2Distrib.md).

The kernel evaluates them through the decomposition described on
[`distrib_deriv3.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Gamma2Distrib.md),
which keeps the digits at a large shape.

## Arguments

- distrib:

  A `Gamma2Distrib` object, from
  [`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md).

- y:

  A numeric vector of strictly positive observations.

- theta:

  A named list with components `mu` and `sigma2`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). The
  generic reads it, and this method does not.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  a fixed number of observations, set per kernel, the kernel runs
  sequentially. Defaults to `1L`.

## Value

A named list of three numeric vectors, `mu_mu`, `sigma2_sigma2` and
`mu_sigma2`, in that order, each of length
`max(length(y), length(mu), length(sigma2))`. The three name the
distinct entries of a symmetric \\2 \times 2\\ matrix per observation.

## Notation

\\\ell^{(ij)}\\ is the second derivative of the log-density in
parameters \\i\\ and \\j\\; parenthesized superscripts name derivatives.
\\\psi\\ and \\\psi_1\\ are the digamma and trigamma functions.

## See also

[`distrib_gradient.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gamma2Distrib.md)
for the score,
[`distrib_expected_hessian.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gamma2Distrib.md)
for the expectation of this quantity,
[`distrib_deriv3.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Gamma2Distrib.md)
for the order above, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- gamma2_distrib()
y <- c(1, 3, 5)
th <- list(mu = 3, sigma2 = 2)
h <- distrib_hessian(d, y, th)
h
#> $mu_mu
#> [1] -1.7219317 -0.6233195 -0.1124938
#> 
#> $sigma2_sigma2
#> [1] -0.8468339  0.1250437 -0.2255986
#> 
#> $mu_sigma2
#> [1]  1.15400317  0.00608474 -0.26015370
#> 

# The curvature in the variance is positive at y = 3, so the observed
# information is not positive definite at every observation.
h$sigma2_sigma2
#> [1] -0.8468339  0.1250437 -0.2255986

# It is the second derivative of the log-density, so a central difference
# of the score reproduces it.
eps <- 1e-6
up <- distrib_gradient(d, y, list(mu = 3 + eps, sigma2 = 2))$mu
dn <- distrib_gradient(d, y, list(mu = 3 - eps, sigma2 = 2))$mu
all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-5)
#> [1] TRUE
```
