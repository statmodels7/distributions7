# Poisson-Inverse Gaussian Third Derivatives

Returns the exact third derivatives of the log-mass in \\(\mu,
\sigma)\\, read off columns `d30`, `d21`, `d12` and `d03` of the
compiled fourth-order kernel of
[`pig_hd_block()`](https://statmodels7.github.io/distributions7/reference/pig_hd_block.md).
The kernel computes all four orders in one pass, so this order costs no
more than the score does.

With `expected = TRUE` the value is an expectation instead, and there it
is **not** closed form: the call routes to
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md),
so `approx` and `nsim` are read.

## Arguments

- distrib:

  A `Pig1Distrib` object, from
  [`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md).

- y:

  A numeric vector of counts. With `expected = TRUE` its values are the
  support points the expectation is summed over.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`.

- expected:

  Logical of length 1. When `FALSE`, the default, the observed
  derivatives at `y` are returned from the compiled kernel.

- scale:

  Either `"parameter"`, the default, or `"link"`. The transformation is
  applied in the generic's body.

- approx:

  One of `"integrate"`, `"bartlett"`, `"mc"` or `"opg"`, the strategy
  for the expectation. Read only when `expected = TRUE`; for a discrete
  family `"integrate"` is an exact sum over the support, not a
  quadrature.

- nsim:

  A single positive integer, the Monte Carlo sample size used when
  `approx = "mc"`. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the compiled kernel may
  use. Defaults to `1L`.

## Value

A named list of four numeric vectors: `mu_mu_mu`, `mu_mu_sigma`,
`mu_sigma_sigma` and `sigma_sigma_sigma`.

## See also

[`distrib_hessian.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Pig1Distrib.md)
for the order below,
[`distrib_deriv4.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Pig1Distrib.md)
for the order above,
[`pig_hd_block()`](https://statmodels7.github.io/distributions7/reference/pig_hd_block.md)
for the kernel, and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- pig1_distrib()
y <- 0:6
th <- list(mu = 3, sigma = 0.8)
d3 <- distrib_deriv3(d, y, th)
names(d3)
#> [1] "mu_mu_mu"          "mu_mu_sigma"       "mu_sigma_sigma"   
#> [4] "sigma_sigma_sigma"

# Against a central difference of the analytic Hessian.
eps <- 1e-5
rbind(analytic = d3$mu_mu_sigma,
      numeric = (distrib_hessian(d, y, list(mu = 3, sigma = 0.8 + eps))$mu_mu -
                 distrib_hessian(d, y, list(mu = 3, sigma = 0.8 - eps))$mu_mu) /
                (2 * eps))
#>                 [,1]          [,2]       [,3]       [,4]      [,5]      [,6]
#> analytic -0.01728057 -0.0008797331 0.02767238 0.06346716 0.1012857 0.1384404
#> numeric  -0.01728057 -0.0008797331 0.02767238 0.06346716 0.1012857 0.1384404
#>               [,7]
#> analytic 0.1743171
#> numeric  0.1743171

# The expected version is not closed form and reads the strategy.
distrib_deriv3(d, 0:200, th, expected = TRUE)$mu_mu_mu[1]
#> [1] 0.1016449
```
