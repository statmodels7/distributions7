# Multivariate Student t Third and Fourth Response Derivatives

Compute the third and fourth derivatives of the log-density in the
response, one array per observation. With \\w = \Sigma^{-1}(y-\mu)\\,
\\q = (y-\mu)^\top w\\, \\s = \nu + q\\, \\N = \nu + p\\ and \\\Omega =
\Sigma^{-1}\\, the log-density is \\-\tfrac{N}{2}\log s\\ plus terms
free of \\y\\, and \$\$\ell\_{ijk} =
\frac{2N}{s^2}\left(\Omega\_{ij}w_k + \Omega\_{ik}w_j +
\Omega\_{jk}w_i\right) - \frac{8N}{s^3}\\w_iw_jw_k,\$\$ \$\$\ell\_{ijkl}
= \frac{2N}{s^2}\left(\Omega\_{ij}\Omega\_{kl} +
\Omega\_{ik}\Omega\_{jl} + \Omega\_{il}\Omega\_{jk}\right) -
\frac{8N}{s^3}\sum\_{(xy)(zt)}\Omega\_{xy}w_zw_t +
\frac{48N}{s^4}\\w_iw_jw_kw_l,\$\$ where the middle sum runs over the
six ways of choosing the pair \\(x, y)\\ from \\\\i, j, k, l\\\\, the
other two indices taking the \\w\\s. Each follows from the order below
by \\\partial s/\partial y = 2w\\ and \\\partial w/\partial y =
\Omega\\.

## Arguments

- distrib:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object, from
  [`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md)
  or
  [`mvstudent_t2_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations. A vector of length
  \\p\\ is read as a single observation.

- theta:

  A named list of parameters, each component a single number.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

[`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
a \\p \times p \times p \times n\\ numeric array,
[`distrib_deriv4_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
a \\p \times p \times p \times p \times n\\ one; slice \\i\\ in the last
dimension belongs to row \\i\\ of `y`.

## Details

A multivariate Student t prior over a block of coefficients has a
Hessian that moves with the coefficients, so the derivative of a
marginal criterion in the hyperparameters reads these arrays through
`penalties7::penalty_dhessian_beta()` and
`penalties7::penalty_d2hessian_beta()`.

## See also

[`distrib_hess_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.MvStudentTDistrib.md)
for the order below,
[`distrib_cross3_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross3_y.MvStudentTDistrib.md)
for the derivative in the parameters, and
[`distrib_deriv3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.md)
for the generics.

## Examples

``` r
d <- mvstudent_t1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
set.seed(1)
y <- distrib_rng(d, 3, theta)
d3 <- distrib_deriv3_y(d, y, theta)
dim(d3)
#> [1] 2 2 2 3

# Against a difference of the response Hessian along the first coordinate.
h <- 1e-5
yp <- y; yp[, 1] <- yp[, 1] + h
ym <- y; ym[, 1] <- ym[, 1] - h
max(abs(d3[, , 1, ] -
        (distrib_hess_y(d, yp, theta) - distrib_hess_y(d, ym, theta)) / (2 * h)))
#> [1] 1.595502e-11
```
