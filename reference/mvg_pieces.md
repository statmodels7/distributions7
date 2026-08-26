# The Pieces a Multivariate Gaussian Evaluates From

Assembles, once per call, the mean, the covariance, the inverse
covariance and the log-determinant from a flat parameter vector,
together with the matrix parametrization's derivative arrays when they
are asked for. Every method of
[MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
calls this first and works from the result, so a parameter vector is
unpacked and a matrix factorized once per call instead of once per
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

  Logical of length 1. When `TRUE` the first derivative arrays of the
  covariance are computed and returned as `a`. Defaults to `FALSE`.

- derivs2:

  Logical of length 1. When `TRUE` the second derivative arrays are
  computed as well and returned as `a2`; the first derivatives are
  computed with them, so `derivs` need not also be set. Defaults to
  `FALSE`.

## Value

A named list with `mu` (numeric of length \\p\\), `sigma` and
`sigma_inv` (\\p \times p\\ matrices), `logdet` (the log-determinant of
\\\Sigma\\, a single number), `eta` (the matrix parametrization's free
vector), `p` and `s` (the parametrization itself), plus `a` and `a2`
when asked for: lists of \\p \times p\\ matrices, `a` one per free value
and `a2` one per unordered pair in
[`parameters7::param_tuple_indices()`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.html)
order.

## Details

The arithmetic downstream is written in the covariance whichever side
the parametrization carries, so a precision parametrization is inverted
here. Three things then change relative to a covariance parametrization:
the log-determinant takes the opposite sign, the first derivatives are
carried across by \\\partial\Sigma/\partial\eta_k = -\Sigma A_k
\Sigma\\, and the second derivatives by differentiating that expression
once more, \$\$\frac{\partial^2\Sigma}{\partial\eta_k\partial\eta_l} =
\Sigma\left(A_l\Sigma A_k + A_k \Sigma A_l -
\frac{\partial^2\Omega}{\partial\eta_k\partial\eta_l}\right)\Sigma,\$\$
with \\A_k = \partial\Omega/\partial\eta_k\\. Those conversions are why
the precision form is not the cheaper one here, contrary to what a
reader might expect from the density alone.

## Notation

\\\Sigma\\ is the covariance, \\\Omega = \Sigma^{-1}\\ the precision,
\\\eta\\ the free vector of the matrix parametrization and \\A_k =
\partial M/\partial\eta_k\\ the derivative of whichever matrix the
parametrization carries.

## See also

[`mvgaussian_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian_distrib.md)
for the family and
[`distrib_gradient.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MvGaussianDistrib.md)
for the first consumer.

## Examples

``` r
d <- mvgaussian_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
pc <- distributions7:::mvg_pieces(d, theta, derivs = TRUE)
names(pc)
#> [1] "mu"        "sigma"     "sigma_inv" "logdet"    "eta"       "p"        
#> [7] "s"         "a"        

# sigma_inv really is the inverse, and logdet its log-determinant.
all.equal(pc$sigma %*% pc$sigma_inv, diag(2), check.attributes = FALSE)
#> [1] TRUE
all.equal(pc$logdet, log(det(pc$sigma)))
#> [1] TRUE

# One derivative array per free value of the parametrization.
length(pc$a)
#> [1] 3
round(pc$a[[3]], 4)
#>        [,1]   [,2]
#> [1,] 0.0000 1.1052
#> [2,] 1.1052 0.8000

# The same pieces from the precision side describe the same covariance.
o <- mvgaussian_distrib(2, omega = parameters7::log_cholesky(2))
th_o <- list(mu1 = 0.5, mu2 = -0.3, omega_log_L1 = 0.1,
             omega_log_L2 = -0.2, omega_L2.1 = 0.4)
po <- distributions7:::mvg_pieces(o, th_o)
all.equal(po$logdet, -parameters7::param_logdet(o@param, unlist(th_o)[3:5]))
#> [1] TRUE
```
