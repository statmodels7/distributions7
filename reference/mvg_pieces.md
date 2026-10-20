# The Pieces From Which a Multivariate Gaussian Evaluates

Assembles, once per call, the mean, the covariance, the inverse
covariance and the log-determinant from a flat parameter vector,
together with the derivative arrays of \\\Sigma\\ and of \\\Sigma^{-1}\\
when they are requested. The density, generator, derivative and moment
methods of
[MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
call this first and work from the result, so a parameter vector is
unpacked and a matrix factorized once per call and not once per
component.

## Usage

``` r
mvg_pieces(distrib, theta, derivs = FALSE, derivs2 = FALSE)
```

## Arguments

- distrib:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object.

- theta:

  A named list of parameters, already aligned by the generic, or any
  list whose components are in `distrib@params` order.

- derivs:

  Logical of length 1. When `TRUE` the first derivative arrays are
  computed. Defaults to `FALSE`.

- derivs2:

  Logical of length 1. When `TRUE` the second derivative arrays are
  computed, and the first ones with them, so `derivs` need not also be
  set. Defaults to `FALSE`.

## Value

A named list with `mu` (numeric of length \\p\\), `eta` (the matrix
parametrization's free vector), `p`, `s` (the parametrization itself)
and the elements of
[`mv_matrix_pieces()`](https://statmodels7.github.io/distributions7/reference/mv_matrix_pieces.md):
`sigma` and `sigma_inv` (\\p \times p\\ matrices), `logdet` (the
log-determinant of \\\Sigma\\) and `inverted`; when derivatives are
requested, the list also holds `a` and `ai` (the first derivative arrays
of \\\Sigma\\ and of \\\Sigma^{-1}\\, one per free value) and `dlogdet`;
when `derivs2` is set, it also holds `a2`, `ai2` (one array per
unordered pair in
[`parameters7::param_tuple_indices()`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.html)
order) and `d2logdet`.

## Details

The matrix quantities come from
[`mv_matrix_pieces()`](https://statmodels7.github.io/distributions7/reference/mv_matrix_pieces.md),
which the Student t reads as well. Whichever side the parametrization
carries, the result holds \\\Sigma\\, \\\Sigma^{-1}\\, the
log-determinant of \\\Sigma\\ and the derivatives of both matrices.

## Notation

\\\Sigma\\ is the covariance, \\\Omega = \Sigma^{-1}\\ the precision and
\\\eta\\ the free vector of the matrix parametrization.

## See also

[`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md)
for the family,
[`mv_matrix_pieces()`](https://statmodels7.github.io/distributions7/reference/mv_matrix_pieces.md)
for the matrix quantities, and
[`distrib_gradient.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MvGaussianDistrib.md)
for a consumer.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
pc <- distributions7:::mvg_pieces(d, theta, derivs = TRUE)
names(pc)
#>  [1] "mu"        "eta"       "p"         "s"         "sigma"     "sigma_inv"
#>  [7] "logdet"    "inverted"  "ai"        "a"         "dlogdet"  

# sigma_inv is the inverse, and logdet its log-determinant.
all.equal(pc$sigma %*% pc$sigma_inv, diag(2), check.attributes = FALSE)
#> [1] TRUE
all.equal(pc$logdet, log(det(pc$sigma)))
#> [1] TRUE

# There is one derivative array per free value of the parametrization.
length(pc$a)
#> [1] 3
round(pc$a[[3]], 4)
#>        [,1]   [,2]
#> [1,] 0.0000 1.1052
#> [2,] 1.1052 0.8000

# The same pieces from the precision side describe the same covariance.
o <- mvgaussian2_distrib(2, parameters7::log_cholesky(2))
th_o <- list(mu1 = 0.5, mu2 = -0.3, omega_log_L1 = 0.1,
             omega_log_L2 = -0.2, omega_L2.1 = 0.4)
po <- distributions7:::mvg_pieces(o, th_o)
all.equal(po$logdet, -parameters7::param_logdet(o@param, unlist(th_o)[3:5]))
#> [1] TRUE
```
