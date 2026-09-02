# Multivariate Student t Third Derivative in One Response and Two Parameters

Computes \\\partial^3\ell/\partial
y\\\partial\theta_a\partial\theta_b\\, one \\n \times p\\ matrix per
unordered pair of parameters. The response gradient is \\\ell^{(y)} =
-c\\w\\, so a second derivative in the parameters is the ordinary
product rule \$\$\partial\_{ab}\ell^{(y)} = -c\_{ab}w - c_aw_b -
c_bw_a - c\\w\_{ab},\$\$ with every piece from
[`mvt_dpieces()`](https://statmodels7.github.io/distributions7/reference/mvt_dpieces.md).
No pair vanishes here, where a gaussian's pair of location parameters is
exactly zero: its response gradient is linear in the location and this
one is not.

## Arguments

- distrib:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object, from
  [`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations. A vector of length
  \\p\\ is read as a single observation.

- theta:

  A named list of parameters, each component a single number.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch. The two differ here.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of \\n \times p\\ numeric matrices, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## Notation

\\\nu\\ is the degrees of freedom, \\c\\ the weight, \\w =
\Sigma^{-1}(y-\mu)\\, and a subscript \\a\\ or \\ab\\ denotes a
derivative in the parameters named.

## See also

[`distrib_cross_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.MvStudentTDistrib.md)
for the order below,
[`distrib_hess_y_hess.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.MvStudentTDistrib.md)
for the sibling with two response indices,
[`mvt_dpieces()`](https://statmodels7.github.io/distributions7/reference/mvt_dpieces.md)
for the assembly, and
[`distrib_grad_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.md)
for the generic.

## Examples

``` r
d <- mvstudent_t1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
set.seed(1)
y <- distrib_rng(d, 4, theta)

gh <- distrib_grad_y_hess(d, y, theta)
dim(gh$mu1_mu2)
#> [1] 4 2

# A gaussian's two-location pair is exactly zero and this one is not.
gh$mu1_mu2
#>            [,1]        [,2]
#> [1,]  0.3831362 -0.51490494
#> [2,] -0.5298276  0.55887216
#> [3,]  0.3837536 -0.62965640
#> [4,] -0.1588397  0.01859468
distrib_grad_y_hess(mvgaussian1_distrib(2), y, theta[1:5])$mu1_mu2
#>      [,1] [,2]
#> [1,]    0    0
#> [2,]    0    0
#> [3,]    0    0
#> [4,]    0    0

# Against a second difference of the response gradient. The reference
# amplifies rounding by h^-2, so 1e-6 on a quantity of order 1 is its own
# accuracy rather than a disagreement.
h <- 1e-5
f <- function(a, b) {
  t2 <- theta; t2$mu1 <- t2$mu1 + a; t2$mu2 <- t2$mu2 + b
  distrib_grad_y(d, y, t2)
}
c(gap = max(abs(gh$mu1_mu2 -
                (f(h, h) - f(h, -h) - f(-h, h) + f(-h, -h)) / (4 * h * h))),
  scale = max(abs(gh$mu1_mu2)))
#>          gap        scale 
#> 1.161690e-06 6.296564e-01 
```
