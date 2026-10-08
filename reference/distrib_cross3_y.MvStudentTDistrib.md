# Multivariate Student t Third Response Derivative in a Parameter

Computes \\\partial^4\ell/\partial y_i\partial y_j\partial y_k\\
\partial\theta_a\\, one \\p \times p \times p \times n\\ array per
parameter. With the notation of
[`distrib_deriv3_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.MvStudentTDistrib.md),
write \\T\_{ijk} = \Omega\_{ij}w_k + \Omega\_{ik}w_j +
\Omega\_{jk}w_i\\. The derivative of \\\ell\_{ijk} = 2Ns^{-2}T\_{ijk} -
8Ns^{-3}w_iw_jw_k\\ in \\\theta_a\\ is \$\$\left(\frac{2N_a}{s^2} -
\frac{4Ns_a}{s^3}\right)T\_{ijk} + \frac{2N}{s^2}\\\partial_aT\_{ijk} -
\left(\frac{8N_a}{s^3} - \frac{24Ns_a}{s^4}\right)w_iw_jw_k -
\frac{8N}{s^3}\\\partial_a(w_iw_jw_k),\$\$ with \\N_a\\, \\s_a\\,
\\\partial_aw\\ and \\\partial_a\Omega\\ from
[`mvt_dpieces()`](https://statmodels7.github.io/distributions7/reference/mvt_dpieces.md):
\\N_a\\ is one for \\\nu\\ and zero otherwise.

## Arguments

- distrib:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object.

- y:

  An \\n \times p\\ numeric matrix of observations.

- theta:

  A named list of parameters, each component a single number.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of \\p \times p \times p \times n\\ numeric arrays, one per
parameter, in `distrib@params` order.

## See also

[`distrib_deriv3_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_y.MvStudentTDistrib.md),
whose derivative in the parameters this is, and
[`distrib_cross3_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross3_y.md)
for the generic.

## Examples

``` r
d <- mvstudent_t1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
set.seed(1)
y <- distrib_rng(d, 3, theta)
c3 <- distrib_cross3_y(d, y, theta)

# Against a difference of the third response derivative in nu.
h <- 1e-5
tp <- theta; tp$nu <- 6 + h
tm <- theta; tm$nu <- 6 - h
max(abs(c3$nu - (distrib_deriv3_y(d, y, tp) -
                 distrib_deriv3_y(d, y, tm)) / (2 * h)))
#> [1] 1.882004e-11
```
