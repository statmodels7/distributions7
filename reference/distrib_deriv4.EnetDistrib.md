# Elastic-Net Fourth Derivatives

Computes the fifteen fourth derivatives of the log-density in closed
form, through
[`.enet_chain()`](https://statmodels7.github.io/distributions7/reference/dot-enet_chain.md),
in the notation of
[`distrib_deriv3.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.EnetDistrib.md).
At this order **nothing touches the location**: the data term is
quadratic in \\z\\ and linear in each rate, so every component naming
\\\mu\\ more than twice is exactly zero, and what remains is a
derivative of the normalizing constant carried through the bilinear map.

With `expected = TRUE` the value is an expectation and is not closed
form, as at third order.

## Arguments

- distrib:

  An `EnetDistrib` object, from
  [`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md).

- y:

  A numeric vector of observations. With `expected = TRUE` only its
  length matters.

- theta:

  A named list with components `mu`, `lambda` and `alpha`, each a
  numeric vector of length 1 or of the length of `y`.

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

## Value

A named list of fifteen numeric vectors, from `mu_mu_mu_mu` to
`alpha_alpha_alpha_alpha`.

## See also

[`distrib_deriv3.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.EnetDistrib.md)
for the order below,
[`.enet_chain()`](https://statmodels7.github.io/distributions7/reference/dot-enet_chain.md)
for the assembly, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- enet_distrib()
y <- c(-1.5, -0.3, 0.4, 2.1)
th <- list(mu = 0, lambda = 2, alpha = 0.5)
d4 <- distrib_deriv4(d, y, th)
length(d4)
#> [1] 15

# Nothing at this order touches the location.
d4$mu_mu_mu_mu
#> [1] 0 0 0 0

# Against a central difference of the third order.
eps <- 1e-5
rbind(analytic = d4$mu_mu_lambda_lambda,
      numeric = (distrib_deriv3(d, y, list(mu = 0, lambda = 2 + eps,
                                           alpha = 0.5))$mu_mu_lambda -
                 distrib_deriv3(d, y, list(mu = 0, lambda = 2 - eps,
                                           alpha = 0.5))$mu_mu_lambda) /
                (2 * eps))
#>          [,1] [,2] [,3] [,4]
#> analytic    0    0    0    0
#> numeric     0    0    0    0

# The pure-lambda component against one stencil on the log-density.
ld <- function(v) sum(distrib_pdf(d, y, list(mu = 0, lambda = v,
                                             alpha = 0.5), log = TRUE))
c(ours = sum(d4$lambda_lambda_lambda_lambda),
  stencil = numericals7::fd_derivative(ld, 2, 4L, h = 0.05))
#>       ours    stencil 
#> -0.9319415 -0.9338097 
```
