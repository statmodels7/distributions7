# Orthogonal Poisson-Inverse Gaussian Fourth Derivatives

Returns the fourth derivatives of the log-mass in \\(\mu, \alpha)\\,
written out in closed form, read off the last five columns of the
compiled kernel `pig2_deriv4_cpp`, which computes the fourth order
alone.

At a large count this order loses relative accuracy for the reason given
on
[`distrib_deriv4.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Pig1Distrib.md):
at \\\alpha = 1\\ and \\y = 400\\ the component in \\\alpha\\ alone has
a relative error of about \\10^{-5}\\. The orders below, and this order
at counts under 100, are accurate to \\10^{-8}\\ or better.

With `expected = TRUE` the value is an expectation and is not closed
form, as at third order.

## Arguments

- distrib:

  A `Pig2Distrib` object, from
  [`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md).

- y:

  A numeric vector of counts. With `expected = TRUE` only its length is
  used, and the expectation is summed over the whole support.

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

A named list of five numeric vectors: `mu_mu_mu_mu`, `mu_mu_mu_alpha`,
`mu_mu_alpha_alpha`, `mu_alpha_alpha_alpha` and
`alpha_alpha_alpha_alpha`.

## See also

[`distrib_deriv3.Pig2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Pig2Distrib.md)
for the order below, `pig2_deriv4_cpp` for the kernel, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- pig2_distrib()
y <- 0:6
al <- sqrt(1 / 0.8^2 + 2 * 3 / 0.8)
th <- list(mu = 3, alpha = al)
d4 <- distrib_deriv4(d, y, th)
names(d4)
#> [1] "mu_mu_mu_mu"             "mu_mu_mu_alpha"         
#> [3] "mu_mu_alpha_alpha"       "mu_alpha_alpha_alpha"   
#> [5] "alpha_alpha_alpha_alpha"

# Against a central difference of the third order.
eps <- 1e-5
rbind(analytic = d4$mu_mu_mu_alpha,
      numeric = (distrib_deriv3(d, y,
                   list(mu = 3, alpha = al + eps))$mu_mu_mu -
                 distrib_deriv3(d, y,
                   list(mu = 3, alpha = al - eps))$mu_mu_mu) / (2 * eps))
#>                [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
#> analytic 0.00993897 0.01965257 0.02936616 0.03907976 0.04879335 0.05850695
#> numeric  0.00993897 0.01965257 0.02936616 0.03907976 0.04879335 0.05850695
#>                [,7]
#> analytic 0.06822055
#> numeric  0.06822055
```
