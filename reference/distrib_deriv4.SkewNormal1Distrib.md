# Skew Normal Fourth Derivatives

Computes the fifteen fourth derivatives of the log-density in closed
form, in a compiled kernel, in the notation of
[`distrib_deriv3.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.SkewNormal1Distrib.md).
The recursion \\R' = -R(t+R)\\ keeps every one of them a polynomial in
\\t\\ and the inverse Mills ratio, so the fourth order needs nothing the
third did not.

With `expected = TRUE` the value is a numerical expectation, for the
reason it is at third order: the integrals are the ones that block the
expected information.

## Arguments

- distrib:

  A `SkewNormal1Distrib` object, from
  [`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md).

- y:

  A numeric vector of observations. With `expected = TRUE` only its
  length matters.

- theta:

  A named list with components `mu`, `sigma` and `alpha`, each a numeric
  vector of length 1 or of the length of `y`.

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
  use. Defaults to `1L`. The result does not depend on the count.

## Value

A named list of fifteen numeric vectors, one per distinct fourth-order
component, named as
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
names them, from `mu_mu_mu_mu` to `alpha_alpha_alpha_alpha`.

## Notation

\\z = (y-\mu)/\sigma\\, \\t = \alpha z\\ and \\R\\ the inverse Mills
ratio.

## See also

[`distrib_deriv3.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.SkewNormal1Distrib.md)
for the order below,
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic, and
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
for the family.

## Examples

``` r
d <- skewnormal1_distrib()
y <- c(-1.5, -0.3, 0.4, 2.1)
th <- list(mu = 0, sigma = 1, alpha = 3)
d4 <- distrib_deriv4(d, y, th)
length(d4)
#> [1] 15

# Against a central difference of the third order.
eps <- 1e-4
rbind(analytic = d4$mu_mu_alpha_alpha,
      numeric = (distrib_deriv3(d, y, list(mu = 0, sigma = 1, alpha = 3 + eps))$mu_mu_alpha -
                 distrib_deriv3(d, y, list(mu = 0, sigma = 1, alpha = 3 - eps))$mu_mu_alpha) /
                (2 * eps))
#>               [,1]      [,2]      [,3]          [,4]
#> analytic -2.030997 -1.959839 0.7048172 -7.884540e-06
#> numeric  -2.030997 -1.959839 0.7048172 -7.884542e-06

# The count is bit for bit independent of the thread count.
identical(distrib_deriv4(d, y, th, threads = 1L),
          distrib_deriv4(d, y, th, threads = 2L))
#> [1] TRUE
```
