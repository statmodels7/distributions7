# Orthogonal Poisson-Inverse Gaussian Observed Hessian

Returns the exact second derivatives of the log-mass in \\(\mu,
\alpha)\\, read off columns `d20`, `d02` and `d11` of the compiled
kernel `pig2_hessian_cpp`.

The observed mixed entry is not zero in general, but its expectation is
zero. At \\\mu = 3\\ and \\\alpha = 3.0104\\ the entry is proportional
to \\y - \mu\\, and its expectation is zero to rounding error.

## Arguments

- distrib:

  A `Pig2Distrib` object, from
  [`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md).

- y:

  A numeric vector of counts. A value off the support, or a missing
  value, gives `NaN`.

- theta:

  A named list with components `mu` and `alpha`, each a numeric vector
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
order: `mu_mu`, `alpha_alpha`, `mu_alpha`.

## See also

[`distrib_gradient.Pig2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Pig2Distrib.md)
for the order below,
[`distrib_deriv3.Pig2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Pig2Distrib.md)
for the order above, `pig2_hessian_cpp` for the kernel, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- pig2_distrib()
y <- 0:6
al <- sqrt(1 / 0.8^2 + 2 * 3 / 0.8)
th <- list(mu = 3, alpha = al)
h <- distrib_hessian(d, y, th)
names(h)
#> [1] "mu_mu"       "alpha_alpha" "mu_alpha"   

# The mixed entry is compared with a central difference of the score.
eps <- 1e-5
rbind(analytic = h$mu_alpha,
      numeric = (distrib_gradient(d, y, list(mu = 3, alpha = al + eps))$mu -
                 distrib_gradient(d, y, list(mu = 3, alpha = al - eps))$mu) /
                (2 * eps))
#>                [,1]       [,2]        [,3]          [,4]       [,5]      [,6]
#> analytic -0.1176464 -0.0784309 -0.03921545 -2.220446e-16 0.03921545 0.0784309
#> numeric  -0.1176464 -0.0784309 -0.03921545  0.000000e+00 0.03921545 0.0784309
#>               [,7]
#> analytic 0.1176464
#> numeric  0.1176464

# The mixed entry is not zero observation by observation; its expectation
# is.
c(observed = h$mu_alpha[1],
  expected = distrib_expected_hessian(d, 0, th)$mu_alpha)
#>      observed      expected 
#> -1.176464e-01 -1.277618e-17 
```
