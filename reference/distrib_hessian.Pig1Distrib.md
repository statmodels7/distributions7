# Poisson-Inverse Gaussian Observed Hessian

Returns the exact second derivatives of the log-mass in \\(\mu,
\sigma)\\, read off columns `d20`, `d02` and `d11` of the compiled
fourth-order kernel of
[`pig_hd_block()`](https://statmodels7.github.io/distributions7/reference/pig_hd_block.md).

This is the **observed** curvature at the data. The expected information
has no closed form for this family and comes from
[`expected_derivative_methods()`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md),
whose default here is the exact sum over the support; see
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md).

## Arguments

- distrib:

  A `Pig1Distrib` object, from
  [`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md).

- y:

  A numeric vector of counts. A value off the support gives `NaN`.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. Both must be strictly positive.

- scale:

  Either `"parameter"`, the default, or `"link"`. The transformation is
  applied in the generic's body.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the compiled kernel may
  use. Defaults to `1L`.

## Value

A named list of three numeric vectors in
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)'s
order: `mu_mu`, `sigma_sigma`, `mu_sigma`.

## See also

[`distrib_gradient.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Pig1Distrib.md)
for the order below,
[`distrib_deriv3.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Pig1Distrib.md)
for the order above,
[`pig_hd_block()`](https://statmodels7.github.io/distributions7/reference/pig_hd_block.md)
for the kernel, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- pig1_distrib()
y <- 0:6
th <- list(mu = 3, sigma = 0.8)
h <- distrib_hessian(d, y, th)
names(h)
#> [1] "mu_mu"       "sigma_sigma" "mu_sigma"   

# Against a central difference of the score.
eps <- 1e-5
rbind(analytic = h$mu_sigma,
      numeric = (distrib_gradient(d, y, list(mu = 3, sigma = 0.8 + eps))$mu -
                 distrib_gradient(d, y, list(mu = 3, sigma = 0.8 - eps))$mu) /
                (2 * eps))
#>               [,1]      [,2]      [,3]       [,4]        [,5]      [,6]
#> analytic 0.2147728 0.1850463 0.1289895 0.05496181 -0.02676516 -0.109158
#> numeric  0.2147728 0.1850463 0.1289895 0.05496181 -0.02676516 -0.109158
#>                [,7]
#> analytic -0.1890535
#> numeric  -0.1890535

# The curvature in the mean is not of one sign: a count of zero and a
# count in the tail pull it in opposite directions.
distrib_hessian(d, c(0, 3, 30), th)$mu_mu
#> [1]  0.05727274 -0.12548872 -1.08836956
```
