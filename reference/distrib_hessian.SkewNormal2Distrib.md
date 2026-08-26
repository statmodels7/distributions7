# Skew Normal Observed Hessian in the Centered Parametrization

Computes the six second derivatives of the log-density in the centered
parameters, by the second-order chain rule through
[`sn_cp_to_dp()`](https://statmodels7.github.io/distributions7/reference/sn_cp_to_dp.md):
\$\$\dfrac{\partial^2 \ell}{\partial \psi_i \partial \psi_j} =
\sum\_{k,l} \dfrac{\partial^2 \ell}{\partial\theta_k\partial\theta_l}
\dfrac{\partial\theta_k}{\partial\psi_i}
\dfrac{\partial\theta_l}{\partial\psi_j} + \sum\_{k} \dfrac{\partial
\ell}{\partial\theta_k}
\dfrac{\partial^2\theta_k}{\partial\psi_i\partial\psi_j}.\$\$ Both terms
come from
[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md),
with the map's partial derivatives supplied by
[`md_skewnormal2()`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md)
as a written-out table.

The **observed** curvature in \\\gamma_1\\ diverges as the skewness goes
to zero, at the rate \\\gamma_1^{-2/3}\\ the cube root sets. Measured at
\\y = 0.5\\, \\\mu = 0\\, \\\sigma = 1\\, it is 0.29, 11.5, 253 and 5451
at \\\gamma_1 = 10^{-2}, 10^{-4}, 10^{-6}, 10^{-8}\\, a factor of 4.642
per decade against \\10^{2/3} = 4.6416\\. The **expected** curvature
does not: see
[`distrib_expected_hessian.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.SkewNormal2Distrib.md).

## Arguments

- distrib:

  A `SkewNormal2Distrib` object, from
  [`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma` and `gamma1`. The skewness
  must not be exactly zero.

- scale:

  Either `"parameter"`, the default, or `"link"`. The transformation is
  applied in the generic's body.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of six numeric vectors, in
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)'s
order: `mu_mu`, `sigma_sigma`, `gamma1_gamma1`, `mu_sigma`, `mu_gamma1`,
`sigma_gamma1`.

## Errors

Signals an error when any element of `gamma1` is exactly zero.

## See also

[`distrib_gradient.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewNormal2Distrib.md)
for the order below,
[`distrib_expected_hessian.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.SkewNormal2Distrib.md)
for the expectation,
[`sn2_chain()`](https://statmodels7.github.io/distributions7/reference/sn2_chain.md)
for the partition sum, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- skewnormal2_distrib()
y <- c(-1, 0.3, 1.7)
th <- list(mu = 0, sigma = 1, gamma1 = 0.5)
h <- distrib_hessian(d, y, th)
names(h)
#> [1] "mu_mu"         "sigma_sigma"   "gamma1_gamma1" "mu_sigma"     
#> [5] "mu_gamma1"     "sigma_gamma1" 

# Against numerical differentiation of the log-density.
f <- function(p) sum(distrib_pdf(d, y, as.list(setNames(p, names(th))),
                                 log = TRUE))
H <- numDeriv::hessian(f, unlist(th))
rbind(analytic = c(sum(h$mu_gamma1), sum(h$gamma1_gamma1)),
      numeric = c(H[1, 3], H[3, 3]))
#>                [,1]      [,2]
#> analytic -0.3985591 0.5335576
#> numeric  -0.3985591 0.5335576

# The observed curvature in the skewness diverges as gamma1^(-2/3). The
# ratio over two decades converges to 10^(2/3) = 4.6416.
cv <- vapply(10^-c(2, 4, 6, 8),
             function(v) distrib_hessian(d, 0.5,
                           list(mu = 0, sigma = 1, gamma1 = v))$gamma1_gamma1, 0)
rbind(curvature = cv, per_decade = c(NA, (cv[-1] / cv[-4])^(1 / 2)))
#>                [,1]      [,2]       [,3]       [,4]
#> curvature  0.288818 11.513887 252.620070 5451.05209
#> per_decade       NA  6.313917   4.684065    4.64522
```
