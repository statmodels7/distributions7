# Orthogonal Poisson-Inverse Gaussian Score

Returns the exact first derivatives of the log-mass in \\(\mu,
\alpha)\\, read off columns `d10` and `d01` of the compiled kernel of
`pig2_gradient_cpp`. The kernel takes \\\alpha\\ as a variable of its
own, so these are derivatives in the orthogonal coordinates directly and
no chain rule through
[`pig2_sigma()`](https://statmodels7.github.io/distributions7/reference/pig2_sigma.md)
is composed.

The property that gives this parametrization its name lives one order
up: the **expected** value of
\\\partial^2\ell/\partial\mu\partial\alpha\\ is zero, so the two scores
are uncorrelated and their maximum likelihood estimates asymptotically
independent.

## Arguments

- distrib:

  A `Pig2Distrib` object, from
  [`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md).

- y:

  A numeric vector of counts. A value off the support gives `NaN`.

- theta:

  A named list with components `mu` and `alpha`, each a numeric vector
  of length 1 or of the length of `y`. Both must be strictly positive.

- scale:

  Either `"parameter"`, the default, or `"link"`. The transformation to
  the link scale is applied in the generic's body, so this method always
  returns the parameter scale.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the compiled kernel may
  use. Defaults to `1L`.

## Value

A named list of two numeric vectors, `mu` and `alpha`, each of the
length of the recycled inputs.

## See also

[`distrib_hessian.Pig2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Pig2Distrib.md)
for the second derivatives,
[`distrib_gradient.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Pig1Distrib.md)
for the same quantity in mean and dispersion, `pig2_gradient_cpp` for
the kernel, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- pig2_distrib()
y <- 0:6
al <- sqrt(1 / 0.8^2 + 2 * 3 / 0.8)
th <- list(mu = 3, alpha = al)
g <- distrib_gradient(d, y, th)

# Against numerical differentiation of the log-likelihood.
f <- function(p) sum(distrib_pdf(d, y, list(mu = p[1], alpha = p[2]),
                                 log = TRUE))
rbind(analytic = vapply(g, sum, 0),
      numeric = numDeriv::grad(f, c(3, al)))
#>                     mu     alpha
#> analytic -6.661338e-16 0.4698132
#> numeric   5.625652e-12 0.4698132

# The two scores are uncorrelated under the model, which is what
# orthogonality means and what pig1 does not have.
c(pig2 = sum(distrib_expected_hessian(d, 0:200, th,
                                      approx = "bartlett")$mu_alpha),
  pig1 = sum(distrib_expected_hessian(pig1_distrib(), 0:200,
                                      list(mu = 3, sigma = 0.8),
                                      approx = "bartlett")$mu_sigma))
#>          pig2          pig1 
#> -1.456596e-14  7.392208e+00 
```
