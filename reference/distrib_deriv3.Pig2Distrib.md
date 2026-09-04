# Orthogonal Poisson-Inverse Gaussian Third Derivatives

Returns the exact third derivatives of the log-mass in \\(\mu,
\alpha)\\, read off columns `d30`, `d21`, `d12` and `d03` of the
compiled kernel `pig2_deriv3_cpp`. All four orders come out of one pass
of the kernel, so this order costs what the score does.

With `expected = TRUE` the value is an expectation instead, and there it
is **not** closed form: the call routes to
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md),
so `approx` and `nsim` are read.

## Arguments

- distrib:

  A `Pig2Distrib` object, from
  [`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md).

- y:

  A numeric vector of counts. With `expected = TRUE` its values are the
  support points the expectation is summed over.

- theta:

  A named list with components `mu` and `alpha`, each a numeric vector
  of length 1 or of the length of `y`.

- expected:

  Logical of length 1. When `FALSE`, the default, the observed
  derivatives at `y` are returned from the compiled kernel.

- scale:

  Either `"parameter"`, the default, or `"link"`. The transformation is
  applied in the generic's body.

- approx:

  One of `"integrate"`, `"bartlett"`, `"mc"` or `"opg"`, read only when
  `expected = TRUE`; for a discrete family `"integrate"` is an exact sum
  over the support.

- nsim:

  A single positive integer, the Monte Carlo sample size used when
  `approx = "mc"`. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the compiled kernel may
  use. Defaults to `1L`.

## Value

A named list of four numeric vectors: `mu_mu_mu`, `mu_mu_alpha`,
`mu_alpha_alpha` and `alpha_alpha_alpha`.

## See also

[`distrib_hessian.Pig2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Pig2Distrib.md)
for the order below,
[`distrib_deriv4.Pig2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Pig2Distrib.md)
for the order above, `pig2_deriv3_cpp` for the kernel, and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- pig2_distrib()
y <- 0:6
al <- sqrt(1 / 0.8^2 + 2 * 3 / 0.8)
th <- list(mu = 3, alpha = al)
d3 <- distrib_deriv3(d, y, th)
names(d3)
#> [1] "mu_mu_mu"          "mu_mu_alpha"       "mu_alpha_alpha"   
#> [4] "alpha_alpha_alpha"

# Against a central difference of the analytic Hessian.
eps <- 1e-5
rbind(analytic = d3$mu_mu_alpha,
      numeric = (distrib_hessian(d, y, list(mu = 3, alpha = al + eps))$mu_mu -
                 distrib_hessian(d, y, list(mu = 3, alpha = al - eps))$mu_mu) /
                (2 * eps))
#>                [,1]          [,2]        [,3]        [,4]        [,5]
#> analytic 0.01940419 -0.0001356936 -0.01967557 -0.03921545 -0.05875533
#> numeric  0.01940419 -0.0001356936 -0.01967557 -0.03921545 -0.05875533
#>                 [,6]        [,7]
#> analytic -0.07829521 -0.09783509
#> numeric  -0.07829521 -0.09783509
```
