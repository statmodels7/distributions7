# The Pieces Every Mixed Response Derivative of the t Is Written In

Assembles, once per call, the first and second derivatives in the
parameters of the four quantities every mixed response derivative of
this family is built from: \\s = \nu + q\\, the whitened residual \\w\\,
the inverse scale matrix \\\Sigma^{-1}\\, and the two scalars \\c =
(\nu+p)/s\\ and \\d = (\nu+p)/s^2\\ that follow from \\s\\. The response
derivatives are \\\ell^{(y)} = -c\\w\\ and \\\ell^{(yy)} =
-c\\\Sigma^{-1} + 2d\\ww^\top\\, so differentiating either in the
parameters is a product rule over these pieces and nothing else.

## Usage

``` r
mvt_dpieces(distrib, y, theta)
```

## Arguments

- distrib:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object, from
  [`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md).

- y:

  An \\n \times p\\ numeric matrix of observations, already coerced by
  [`as_mv_matrix()`](https://statmodels7.github.io/distributions7/reference/as_mv_matrix.md).

- theta:

  A named list of parameters, each component a single number.

## Value

A named list of the quantities above, INDEXED BY PARAMETER POSITION: the
\\p\\ locations, then the matrix parametrization's free values, then
\\\nu\\. Its components are `w`, `si`, `cw`, `dw` and `q` (the values),
`sa`, `wa`, `Sa`, `ca`, `da` (lists of first derivatives, one per
parameter), `npar`, and `pair(a, b)`, a function returning the
second-order pieces `s`, `w`, `S`, `c` and `d` for one pair of
positions.

## First derivatives

Writing \\A_k\\ for the matrix parametrization's first derivatives,
\\u_j\\ for the \\j\\th column of \\\Sigma^{-1}\\ and \\P_k\\ for
\\\Sigma^{-1}A_k\Sigma^{-1}\\, \$\$\partial\_{\mu_j} s = -2w_j, \quad
\partial\_{\eta_k} s = -w^\top A_kw, \quad \partial\_\nu s = 1,\$\$
\$\$\partial\_{\mu_j} w = -u_j, \quad \partial\_{\eta_k} w =
-\Sigma^{-1}A_kw, \quad \partial\_\nu w = 0, \qquad
\partial\_{\eta_k}\Sigma^{-1} = -P_k.\$\$

## Second derivatives

All vanish except \$\$\partial\_{\mu_j\mu_l} s = 2(\Sigma^{-1})\_{jl},
\quad \partial\_{\mu_j\eta_k} s = 2(\Sigma^{-1}A_kw)\_j, \quad
\partial\_{\mu_j\eta_k} w = P_ke_j,\$\$ and, with \\B\_{kl} =
A_k\Sigma^{-1}A_l + A_l\Sigma^{-1}A_k - A\_{kl}\\,
\$\$\partial\_{\eta_k\eta_l} s = w^\top B\_{kl}w, \quad
\partial\_{\eta_k\eta_l} w = \Sigma^{-1}B\_{kl}w, \quad
\partial\_{\eta_k\eta_l}\Sigma^{-1} = \Sigma^{-1}B\_{kl}\Sigma^{-1}.\$\$
One middle matrix \\B\_{kl}\\ serves all three, which keeps the assembly
short. \\c\\ and \\d\\ then follow by the quotient rule, \\\nu+p\\ being
linear in \\\nu\\ and constant in everything else.

## Notation

\\\mu\\ is the location, \\\Sigma\\ the scale matrix, \\\nu\\ the
degrees of freedom, \\\eta\\ the free vector of the matrix
parametrization, \\A_k\\ and \\A\_{kl}\\ its derivative arrays, \\q\\
the squared Mahalanobis distance, \\w = \Sigma^{-1}(y-\mu)\\ and \\e_j\\
the \\j\\th standard basis vector.

## See also

[`distrib_cross_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.MvStudentTDistrib.md),
[`distrib_cross2_y.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.MvStudentTDistrib.md),
[`distrib_grad_y_hess.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y_hess.MvStudentTDistrib.md)
and
[`distrib_hess_y_hess.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y_hess.MvStudentTDistrib.md),
the four consumers, and
[`mvt_weights()`](https://statmodels7.github.io/distributions7/reference/mvt_weights.md)
for the values these differentiate.

## Examples

``` r
d <- mvstudent_t1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
set.seed(1)
y <- distrib_rng(d, 4, theta)
z <- distributions7:::mvt_dpieces(d, y, theta)
names(z)
#>  [1] "p"    "n"    "npar" "si"   "w"    "s"    "cw"   "dw"   "sa"   "wa"  
#> [11] "Sa"   "ca"   "da"   "pair"

# Six parameters, so six first-derivative entries per quantity.
c(npar = z$npar, n_sa = length(z$sa), n_wa = length(z$wa))
#> npar n_sa n_wa 
#>    6    6    6 

# ds/dnu is exactly one, s being nu + q.
z$sa[[6]]
#> [1] 1 1 1 1

# And ds/dmu1 is -2 w1, against a difference of q + nu.
h <- 1e-5
sfun <- function(m1) {
  t2 <- theta; t2$mu1 <- m1
  pc <- distributions7:::mvt_pieces(d, t2)
  distributions7:::mvt_weights(y, pc)$q + t2$nu
}
max(abs(z$sa[[1]] - (sfun(0.5 + h) - sfun(0.5 - h)) / (2 * h)))
#> [1] 6.164669e-11
```
