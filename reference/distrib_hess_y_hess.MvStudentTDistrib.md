# Multivariate Student t Fourth Derivative in Two Responses and Two Parameters

Computes \\\partial^4\ell/\partial y\\\partial y^\top
\partial\theta_a\partial\theta_b\\, one \\p \times p \times n\\ array
per unordered pair of parameters. It is the expansion of \\M =
\ell^{(yy)} = -c\\\Sigma^{-1} + 2d\\ww^\top\\ carried one order further
than
[`distrib_cross2_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.MvStudentTDistrib.md):
\$\$\partial\_{ab} M = -c\_{ab}\Sigma^{-1} - c_a\partial_b\Sigma^{-1} -
c_b\partial_a\Sigma^{-1} - c\\\partial\_{ab}\Sigma^{-1} +
2d\_{ab}\\ww^\top + 2d_a S(w_b, w) + 2d_b S(w_a, w) +
2d\left\\S(w\_{ab}, w) + S(w_a, w_b)\right\\,\$\$ with \\S(u, v) =
uv^\top + vu^\top\\ and every piece from
[`mvt_dpieces()`](https://statmodels7.github.io/distributions7/reference/mvt_dpieces.md).
No pair vanishes, where for a gaussian every pair naming a location
parameter is exactly zero.

## Arguments

- distrib:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object, from
  [`mvstudent_t_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t_distrib.md).

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

A named list of \\p \times p \times n\\ numeric arrays, keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md).
Every slice is symmetric, the derivative of a symmetric matrix being
symmetric.

## Notation

\\\Sigma\\ is the scale matrix, \\c\\ and \\d\\ the two weights, \\w =
\Sigma^{-1}(y-\mu)\\, and a subscript \\a\\ or \\ab\\ denotes a
derivative in the parameters named.

## See also

[`distrib_cross2_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.MvStudentTDistrib.md)
for the same derivative one parameter down,
[`distrib_grad_y_hess.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.MvStudentTDistrib.md)
for the sibling with one response index,
[`mvt_dpieces()`](https://statmodels7.github.io/distributions7/reference/mvt_dpieces.md)
for the assembly, and
[`distrib_hess_y_hess()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.md)
for the generic.

## Examples

``` r
d <- mvstudent_t_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
set.seed(1)
y <- distrib_rng(d, 4, theta)

hh <- distrib_hess_y_hess(d, y, theta)
dim(hh$nu_nu)
#> [1] 2 2 4

# Every slice is symmetric.
all(vapply(seq_len(4), function(i) isSymmetric(hh$nu_nu[, , i]), TRUE))
#> [1] TRUE

# Against a second difference of the response Hessian, with the reference's
# own h^-2 amplification beside it.
h <- 1e-5
f <- function(a, b) {
  t2 <- theta
  t2$sigma_log_L1 <- t2$sigma_log_L1 + a
  t2$sigma_L2.1 <- t2$sigma_L2.1 + b
  distrib_hess_y(d, y, t2)
}
key <- "sigma_log_L1_sigma_L2.1"
c(gap = max(abs(hh[[key]] -
                (f(h, h) - f(h, -h) - f(-h, h) + f(-h, -h)) / (4 * h * h))),
  scale = max(abs(hh[[key]])))
#>          gap        scale 
#> 7.384823e-07 1.163532e+00 
```
