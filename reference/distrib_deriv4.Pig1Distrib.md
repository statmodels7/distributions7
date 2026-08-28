# Poisson-Inverse Gaussian Fourth Derivatives

Returns the exact fourth derivatives of the log-mass in \\(\mu,
\sigma)\\, read off the last five columns of the compiled kernel of
[`pig_hd_block()`](https://statmodels7.github.io/distributions7/reference/pig_hd_block.md).
The kernel is a fourth-order one throughout, so this is the order it was
written for and it costs the same as the score.

With `expected = TRUE` the value is an expectation and is not closed
form, as at third order.

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

  One of `"integrate"`, `"bartlett"`, `"mc"` or `"opg"`, read only when
  `expected = TRUE`.

- nsim:

  A single positive integer, the Monte Carlo sample size used when
  `approx = "mc"`. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the compiled kernel may
  use. Defaults to `1L`.

## Value

A named list of five numeric vectors: `mu_mu_mu_mu`, `mu_mu_mu_sigma`,
`mu_mu_sigma_sigma`, `mu_sigma_sigma_sigma` and
`sigma_sigma_sigma_sigma`.

## See also

[`distrib_deriv3.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Pig1Distrib.md)
for the order below,
[`pig_hd_block()`](https://statmodels7.github.io/distributions7/reference/pig_hd_block.md)
for the kernel, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- pig1_distrib()
y <- 0:6
th <- list(mu = 3, sigma = 0.8)
d4 <- distrib_deriv4(d, y, th)
names(d4)
#> [1] "mu_mu_mu_mu"             "mu_mu_mu_sigma"         
#> [3] "mu_mu_sigma_sigma"       "mu_sigma_sigma_sigma"   
#> [5] "sigma_sigma_sigma_sigma"

# Against a central difference of the third order.
eps <- 1e-5
rbind(analytic = d4$mu_mu_mu_sigma,
      numeric = (distrib_deriv3(d, y, list(mu = 3, sigma = 0.8 + eps))$mu_mu_mu -
                 distrib_deriv3(d, y, list(mu = 3, sigma = 0.8 - eps))$mu_mu_mu) /
                (2 * eps))
#>                 [,1]        [,2]        [,3]        [,4]        [,5]       [,6]
#> analytic 0.002043023 -0.01153008 -0.03446982 -0.06253726 -0.09171889 -0.1203852
#> numeric  0.002043023 -0.01153008 -0.03446982 -0.06253726 -0.09171889 -0.1203852
#>               [,7]
#> analytic -0.148394
#> numeric  -0.148394

# All four orders come from one pass of the kernel, so the fourth costs
# what the first does.
n <- 2e4
set.seed(62)
x <- distrib_rng(d, n, th)
rbind(score = system.time(distrib_gradient(d, x, th))[["elapsed"]],
      fourth = system.time(distrib_deriv4(d, x, th))[["elapsed"]])
#>         [,1]
#> score  0.017
#> fourth 0.016
```
