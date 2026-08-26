# Elastic-Net Third Derivatives

Computes the ten third derivatives of the log-density in closed form,
through
[`.enet_chain()`](https://statmodels7.github.io/distributions7/reference/dot-enet_chain.md).
Written in \\(\mu, a, c)\\ the log-density is quadratic in \\z = y-\mu\\
and linear in each rate, so at this order only the normalizing constant
contributes, apart from \\\partial^{3}\ell/\partial\mu^{2}\partial c =
-1\\. The bilinear map to \\(\lambda, \alpha)\\ then puts the data back
in: `mu_mu_lambda` and `mu_mu_alpha` are not constant across
observations even though their preimages are.

The license for this order is that the same assembly at orders one and
two reproduces the hand-written score and Hessian, measured at exactly 0
and \\1.1\times10^{-16}\\.

With `expected = TRUE` the value is an expectation, and there it is
**not** closed form: the call routes to
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md),
so `approx` and `nsim` are read.

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

A named list of ten numeric vectors, from `mu_mu_mu` to
`alpha_alpha_alpha` as
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
names them.

## Notation

\\a = \lambda\alpha\\, \\c = \lambda(1-\alpha)\\, \\z = y - \mu\\, and
\\G = \mathrm{d}\log M/\mathrm{d}x\\ with \\M\\ the Mills ratio.

## See also

[`distrib_hessian.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.EnetDistrib.md)
for the order below,
[`distrib_deriv4.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.EnetDistrib.md)
for the order above,
[`.enet_chain()`](https://statmodels7.github.io/distributions7/reference/dot-enet_chain.md)
for the assembly, and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- enet_distrib()
y <- c(-1.5, -0.3, 0.4, 2.1)
th <- list(mu = 0, lambda = 2, alpha = 0.5)
d3 <- distrib_deriv3(d, y, th)
names(d3)
#>  [1] "mu_mu_mu"             "mu_mu_lambda"         "mu_mu_alpha"         
#>  [4] "mu_lambda_lambda"     "mu_lambda_alpha"      "mu_alpha_alpha"      
#>  [7] "lambda_lambda_lambda" "lambda_lambda_alpha"  "lambda_alpha_alpha"  
#> [10] "alpha_alpha_alpha"   

# Against a central difference of the analytic Hessian.
eps <- 1e-5
rbind(analytic = d3$mu_mu_lambda,
      numeric = (distrib_hessian(d, y, list(mu = 0, lambda = 2 + eps,
                                            alpha = 0.5))$mu_mu -
                 distrib_hessian(d, y, list(mu = 0, lambda = 2 - eps,
                                            alpha = 0.5))$mu_mu) / (2 * eps))
#>          [,1] [,2] [,3] [,4]
#> analytic -0.5 -0.5 -0.5 -0.5
#> numeric  -0.5 -0.5 -0.5 -0.5

# The log-density is quadratic in the response, so the third derivative in
# the location alone is exactly zero.
d3$mu_mu_mu
#> [1] 0 0 0 0

# The pure-lambda component against one stencil on the log-density, a
# route that shares no algebra with the assembly.
ld <- function(v) sum(distrib_pdf(d, y, list(mu = 0, lambda = v,
                                             alpha = 0.5), log = TRUE))
c(ours = sum(d3$lambda_lambda_lambda),
  stencil = numericals7::fd_derivative(ld, 2, 3L, h = 0.02))
#>      ours   stencil 
#> 0.6431019 0.6432840 
```
