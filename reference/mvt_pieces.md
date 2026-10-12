# The Pieces From Which a Multivariate t Evaluates

Assembles, once per call, the location, the scale matrix, its inverse,
the log-determinant and the degrees of freedom from a flat parameter
vector, together with the derivative arrays of \\\Sigma\\ and of
\\\Sigma^{-1}\\ when they are requested. Every method of
[MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
that evaluates the density or its derivatives calls this first and works
from the result, so a parameter vector is unpacked and a matrix
factorized once per call and not once per component.

## Usage

``` r
mvt_pieces(distrib, theta, derivs = FALSE, derivs2 = FALSE)
```

## Arguments

- distrib:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object.

- theta:

  A named list of parameters, already aligned by the generic, or any
  list whose components are in `distrib@params` order.

- derivs:

  Logical of length 1. When `TRUE` the first derivative arrays are
  computed. Defaults to `FALSE`.

- derivs2:

  Logical of length 1. When `TRUE` the second derivative arrays are
  computed, and the first ones with them. Defaults to `FALSE`.

## Value

A named list with `mu` (numeric of length \\p\\), `eta` (the matrix
parametrization's free vector), `nu` (a single number), `p`, `s` (the
parametrization itself) and the elements of
[`mv_matrix_pieces()`](https://statmodels7.github.io/distributions7/reference/mv_matrix_pieces.md):
`sigma` and `sigma_inv` (\\p \times p\\ matrices), `logdet` (the
log-determinant of \\\Sigma\\) and `inverted`; with derivatives
requested also `a` and `ai` (the first derivative arrays of \\\Sigma\\
and of \\\Sigma^{-1}\\) and `dlogdet`; with `derivs2` also `a2`, `ai2`
and `d2logdet`.

## Details

The matrix quantities come from
[`mv_matrix_pieces()`](https://statmodels7.github.io/distributions7/reference/mv_matrix_pieces.md),
which the gaussian reads as well: where `distrib@inverted` is `TRUE` the
parametrization carries \\\Sigma^{-1}\\, and the helper returns
\\\Sigma\\ and the derivatives of both matrices whichever side it
carries.

## Notation

\\\mu\\ is the location, \\\Sigma\\ the scale matrix, \\\nu\\ the
degrees of freedom and \\\eta\\ the free vector of the matrix
parametrization.

## See also

[`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md)
for the family,
[`mvt_weights()`](https://statmodels7.github.io/distributions7/reference/mvt_weights.md)
for the quantities computed from these, and
[`mvg_pieces()`](https://statmodels7.github.io/distributions7/reference/mvg_pieces.md)
for the gaussian's twin.

## Examples

``` r
d <- mvstudent_t1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
pc <- distributions7:::mvt_pieces(d, theta, derivs = TRUE)
names(pc)
#>  [1] "mu"        "eta"       "nu"        "p"         "s"         "sigma"    
#>  [7] "sigma_inv" "logdet"    "inverted"  "ai"        "a"         "dlogdet"  

# sigma_inv is the inverse of the scale matrix, and logdet its determinant.
all.equal(pc$sigma %*% pc$sigma_inv, diag(2), check.attributes = FALSE)
#> [1] TRUE
all.equal(pc$logdet, log(det(pc$sigma)))
#> [1] TRUE

# nu is carried apart from the matrix, as its own number.
c(nu = pc$nu, n_free = length(pc$eta), n_deriv = length(pc$a))
#>      nu  n_free n_deriv 
#>       6       3       3 
```
