# Elastic-Net Observed Hessian

Computes the six second derivatives of the log-density in closed form.
The data term of the log-density is **linear in the two rates**, so
every rate block is a second derivative of \\\log Z\\ carried through
the bilinear map \\(\lambda,\alpha) \mapsto (a, c)\\, whose own cross
term contributes the \\\lambda\alpha\\ entry. In the location the
curvature is simply \\-c\\.

That \\-c\\ is the observed curvature and **not** the information. It
misses the point mass \\\mathrm{d}\\\mathrm{sgn}(z)/\mathrm{d}z =
2\delta(z)\\ the absolute value carries at the location, exactly as the
Laplace does; see
[`distrib_expected_hessian.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.EnetDistrib.md)
for what the information is instead.

## Arguments

- distrib:

  An `EnetDistrib` object, from
  [`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `lambda` and `alpha`, each a
  numeric vector of length 1 or of the length of `y`.

- scale:

  Either `"parameter"`, the default, or `"link"`. The transformation is
  applied in the generic's body.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of six numeric vectors in
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)'s
order: `mu_mu`, `lambda_lambda`, `alpha_alpha`, `mu_lambda`, `mu_alpha`,
`lambda_alpha`. Two of them, `lambda_lambda` and `alpha_alpha`, carry no
data at all and are their own expectations.

## Notation

\\a = \lambda\alpha\\, \\c = \lambda(1-\alpha)\\, \\z = y - \mu\\, \\x =
a/\sqrt c\\, and \\Z\\ the normalizing constant.

## See also

[`distrib_gradient.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.EnetDistrib.md)
for the order below,
[`distrib_deriv3.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.EnetDistrib.md)
for the order above,
[`distrib_expected_hessian.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.EnetDistrib.md)
for the information, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- enet_distrib()
y <- c(-1.5, -0.3, 0.4, 2.1)
th <- list(mu = 0, lambda = 2, alpha = 0.5)
h <- distrib_hessian(d, y, th)
names(h)
#> [1] "mu_mu"         "lambda_lambda" "alpha_alpha"   "mu_lambda"    
#> [5] "mu_alpha"      "lambda_alpha" 

# The curvature in the location is -c at every observation, the absolute
# value contributing nothing away from the kink.
c(observed = unique(h$mu_mu), minus_c = -(2 * (1 - 0.5)))
#> observed  minus_c 
#>       -1       -1 

# Two entries carry no data, so they equal their own expectations.
e <- distrib_expected_hessian(d, y, th)
rbind(observed = c(h$lambda_lambda[1], h$alpha_alpha[1], h$mu_mu[1]),
      expected = c(e$lambda_lambda[1], e$alpha_alpha[1], e$mu_mu[1]))
#>                [,1]       [,2]      [,3]
#> observed -0.1702646 -0.1159321 -1.000000
#> expected -0.1702646 -0.1159321 -2.525135

# Against a central difference of the score.
eps <- 1e-5
rbind(analytic = h$lambda_alpha,
      numeric = (distrib_gradient(d, y, list(mu = 0, lambda = 2,
                                             alpha = 0.5 + eps))$lambda -
                 distrib_gradient(d, y, list(mu = 0, lambda = 2,
                                             alpha = 0.5 - eps))$lambda) /
                (2 * eps))
#>                [,1]        [,2]        [,3]      [,4]
#> analytic -0.1304718 -0.01047179 -0.07547179 0.3495282
#> numeric  -0.1304718 -0.01047179 -0.07547179 0.3495282
```
