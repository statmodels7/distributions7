# Skew t Fourth Derivatives

Computes the thirty-five fourth derivatives of the log-density with the
compiled kernel described on
[`distrib_gradient.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewTDistrib.md).

## Arguments

- distrib:

  A `SkewTDistrib` object, from
  [`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md).

- y:

  A numeric vector of observations. With `expected = TRUE` only its
  length matters.

- theta:

  A named list with components `mu`, `sigma`, `alpha` and `nu`.

- expected:

  Logical of length 1. When `FALSE`, the default, the observed
  derivatives at `y` are returned.

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

  A single positive integer, how many threads the kernel may use.

## Value

A named list of thirty-five numeric vectors, one per distinct
fourth-order component, from `mu_mu_mu_mu` to `nu_nu_nu_nu`.

## Details

With `expected = TRUE` the whole order is an expectation and comes from
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md);
`approx` and `nsim` are then read.

## See also

[`distrib_deriv3.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.SkewTDistrib.md)
for the order below,
[`distrib_deriv5.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.SkewTDistrib.md)
for the order above, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- skewt_distrib()
y <- c(-1.5, -0.3, 0.4, 2.1)
th <- list(mu = 0, sigma = 1, alpha = 3, nu = 6)
d4 <- distrib_deriv4(d, y, th)
c(components = length(d4), involving_nu = sum(grepl("nu", names(d4))))
#>   components involving_nu 
#>           35           20 

# The components span five orders of magnitude, so a relative comparison
# on the smallest of them says nothing about the largest.
s <- sort(vapply(d4, function(v) sum(abs(v)), 0), decreasing = TRUE)
s[c(1, 2, length(s) - 1, length(s))]
#> sigma_sigma_sigma_sigma          mu_mu_mu_sigma             nu_nu_nu_nu 
#>            1.268203e+02            7.633429e+01            4.734965e-03 
#>          alpha_nu_nu_nu 
#>            2.024170e-03 

# A closed-form-block component against a difference of the third order.
eps <- 1e-5
rbind(analytic = d4$mu_mu_alpha_alpha,
      numeric = (distrib_deriv3(d, y, list(mu = 0, sigma = 1,
                                           alpha = 3 + eps, nu = 6))$mu_mu_alpha -
                 distrib_deriv3(d, y, list(mu = 0, sigma = 1,
                                           alpha = 3 - eps, nu = 6))$mu_mu_alpha) /
                (2 * eps))
#>                [,1]      [,2]      [,3]         [,4]
#> analytic -0.2739461 0.9096393 0.7661879 -0.007254261
#> numeric  -0.2739461 0.9096393 0.7661879 -0.007254261
```
