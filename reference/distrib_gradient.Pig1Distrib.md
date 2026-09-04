# Poisson-Inverse Gaussian Score

Returns the exact first derivatives of the log-mass in \\(\mu,
\sigma)\\, read off columns `d10` and `d01` of the compiled kernel of
`pig1_gradient_cpp`. Nothing is differenced: the kernel writes every
partial out in closed form from the finite Bessel sum.

## Arguments

- distrib:

  A `Pig1Distrib` object, from
  [`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md).

- y:

  A numeric vector of counts. A value off the support gives `NaN`.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. Both must be strictly positive.

- scale:

  Either `"parameter"`, the default, or `"link"`. The transformation to
  the link scale is applied in the generic's body, so this method always
  returns the parameter scale.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the compiled kernel may
  use. Defaults to `1L`. The result does not depend on the count.

## Value

A named list of two numeric vectors, `mu` and `sigma`, each of the
length of the recycled inputs.

## Notation

\\\ell\\ is the log-mass of one observation, \\\mu\\ the mean and
\\\sigma\\ the dispersion.

## See also

`pig1_gradient_cpp` for the kernel and the closed form it evaluates,
[`distrib_hessian.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Pig1Distrib.md)
for the second derivatives, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- pig1_distrib()
y <- 0:6
th <- list(mu = 3, sigma = 0.8)
g <- distrib_gradient(d, y, th)

# Against numerical differentiation of the log-likelihood.
f <- function(p) sum(distrib_pdf(d, y, list(mu = p[1], sigma = p[2]),
                                 log = TRUE))
rbind(analytic = vapply(g, sum, 0),
      numeric = numDeriv::grad(f, c(3, 0.8)))
#>                 mu     sigma
#> analytic 0.1950793 -1.036359
#> numeric  0.1950793 -1.036359

# The mean and the dispersion are not orthogonal here: the mixed entry of
# the expected information is far from zero. pig2_distrib() removes that.
sum(distrib_expected_hessian(d, 0:200, th, approx = "bartlett")$mu_sigma)
#> [1] 7.392208
```
