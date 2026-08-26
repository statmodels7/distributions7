# Skew t Second Response Derivative

Computes \\\partial^2\ell/\partial y^2 = D'/\sigma^2\\, with \\D' = A' +
Q'B^2 + QB'\\ in the notation of
[`skewt_pieces()`](https://statmodels7.github.io/distributions7/reference/skewt_pieces.md).
It is the same expression as \\\partial^2\ell/\partial\mu^2\\: two signs
cancel where one did not at first order.

Unlike the skew normal's, the value **changes sign**. The Student
\\t\\'s curvature is positive wherever \\\|z\| \> \sqrt{\nu}\\, so the
log-density is convex in the far tail. That convexity makes the family's
score redescend, and its estimates resistant to an outlier.

## Arguments

- distrib:

  A `SkewTDistrib` object, from
  [`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma`, `alpha` and `nu`, each a
  numeric vector of length 1 or of the length of `y`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of the length of the recycled inputs. It is not of one
sign.

## Notation

\\z = (y-\mu)/\sigma\\, \\\nu\\ the degrees of freedom, and \\A'\\,
\\B\\, \\B'\\, \\Q\\, \\Q'\\ are as
[`skewt_pieces()`](https://statmodels7.github.io/distributions7/reference/skewt_pieces.md)
defines them.

## See also

[`distrib_grad_y.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.SkewTDistrib.md)
for the first derivative,
[`distrib_hessian.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.SkewTDistrib.md)
for the parameter curvature it shares an expression with,
[`distrib_hess_y.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.SkewNormal1Distrib.md)
for the case that does not change sign, and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
for the generic.

## Examples

``` r
d <- skewt_distrib()
y <- c(-1.5, -0.3, 0.4, 2.1)
th <- list(mu = 0, sigma = 1, alpha = 3, nu = 6)

# It equals the curvature in the location, without a sign change.
all.equal(distrib_hess_y(d, y, th), distrib_hessian(d, y, th)$mu_mu)
#> [1] TRUE

# Against a central difference of the response derivative.
eps <- 1e-5
rbind(analytic = distrib_hess_y(d, y, th),
      numeric = (distrib_grad_y(d, y + eps, th) -
                 distrib_grad_y(d, y - eps, th)) / (2 * eps))
#>              [,1]      [,2]      [,3]       [,4]
#> analytic 1.593377 -4.367025 -3.711409 -0.1050938
#> numeric  1.593377 -4.367025 -3.711409 -0.1050938

# Convex in the far tail, where the skew normal is concave everywhere.
sn <- skewnormal1_distrib()
far <- c(-20, -8, 8, 20)
rbind(skew_t = distrib_hess_y(d, far, th),
      skew_normal = distrib_hess_y(sn, far, list(mu = 0, sigma = 1, alpha = 3)))
#>                    [,1]       [,2]        [,3]        [,4]
#> skew_t       0.01743002  0.1066761  0.08285533  0.01673175
#> skew_normal -9.99750415 -9.9845354 -1.00000000 -1.00000000
```
