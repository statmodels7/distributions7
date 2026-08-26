# Skew Normal Score in the Centered Parametrization

Computes the three first derivatives of the log-density in the centered
parameters, by carrying
[`distrib_gradient.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewNormal1Distrib.md)'s
through the Jacobian of
[`sn_cp_to_dp()`](https://statmodels7.github.io/distributions7/reference/sn_cp_to_dp.md):
\$\$\dfrac{\partial \ell}{\partial \psi_j} = \sum\_{k} \dfrac{\partial
\ell}{\partial \theta_k} \dfrac{\partial \theta_k}{\partial \psi_j},
\qquad \psi = (\mu, \sigma, \gamma_1),\\ \theta = (\xi, \omega,
\alpha).\$\$

The component in \\\gamma_1\\ stays of order one however small
\\\gamma_1\\ is, although the Jacobian itself grows without bound.
Measured at \\y = 0.5\\, \\\mu = 0\\, \\\sigma = 1\\, the score reads
\\-0.2152, -0.2257, -0.2284, -0.2290\\ at \\\gamma_1 = 10^{-2}, 10^{-4},
10^{-6}, 10^{-8}\\ while \\\partial\alpha/\partial\gamma_1\\ reads
\\5.1, 258\\ at the first two. The divergent parts cancel, and that
cancellation is the reason the parametrization exists.

## Arguments

- distrib:

  A `SkewNormal2Distrib` object, from
  [`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma` and `gamma1`. The skewness
  must not be exactly zero; see the error below.

- scale:

  Either `"parameter"`, the default, or `"link"`. The transformation is
  applied in the generic's body, so this method always returns the
  parameter scale.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of three numeric vectors, `mu`, `sigma` and `gamma1`, each
of the length of the recycled inputs.

## Errors

Signals an error when any element of `gamma1` is exactly zero: the map
runs through a cube root and is not differentiable there. The density is
defined at that point, and
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
carries the same family with ordinary derivatives at symmetry.

## Notation

\\\ell\\ is the log-density of one observation, \\\psi\\ the centered
parameters and \\\theta\\ the direct ones.

## See also

[`sn2_chain()`](https://statmodels7.github.io/distributions7/reference/sn2_chain.md)
for the partition sum this calls,
[`distrib_hessian.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.SkewNormal2Distrib.md)
for the next order, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- skewnormal2_distrib()
y <- c(-1, 0.3, 1.7)
th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
g <- distrib_gradient(d, y, th)

# Against numerical differentiation of the log-density itself.
f <- function(p) sum(distrib_pdf(d, y, as.list(setNames(p, names(th))),
                                 log = TRUE))
rbind(analytic = vapply(g, sum, 0),
      numeric = numDeriv::grad(f, unlist(th)))
#>                mu     sigma    gamma1
#> analytic 0.772631 0.4856841 0.2120403
#> numeric  0.772631 0.4856841 0.2120403

# The score in the skewness stays bounded as the map's Jacobian diverges.
vapply(10^-c(2, 4, 6, 8),
       function(v) distrib_gradient(d, 0.5,
                     list(mu = 0, sigma = 1, gamma1 = v))$gamma1, 0)
#> [1] -0.2152070 -0.2256687 -0.2284080 -0.2290031
```
