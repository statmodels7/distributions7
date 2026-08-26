# The Pieces a Multivariate t Evaluates From

Assembles, once per call, the location, the scale matrix, its inverse,
the log-determinant and the degrees of freedom from a flat parameter
vector, together with the matrix parametrization's derivative arrays
when they are asked for. Every method of
[MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
calls this first and works from the result, so a parameter vector is
unpacked and a matrix factorized once per call instead of once per
component.

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

  Logical of length 1. When `TRUE` the first derivative arrays of the
  scale matrix are computed and returned as `a`. Defaults to `FALSE`.

- derivs2:

  Logical of length 1. When `TRUE` the second derivative arrays are
  computed as well and returned as `a2`; the first derivatives are not
  computed with them, so a caller needing both sets `derivs` too.
  Defaults to `FALSE`.

## Value

A named list with `mu` (numeric of length \\p\\), `eta` (the matrix
parametrization's free vector), `nu` (a single number), `p`, `s` (the
parametrization itself), `sigma` and `sigma_inv` (\\p \times p\\
matrices) and `logdet` (the log-determinant of \\\Sigma\\), plus `a` and
`a2` when asked for: lists of \\p \times p\\ matrices, `a` one per free
value and `a2` one per unordered pair.

## Details

The scale matrix is always parametrized as itself, so there is no
inversion branch here and no sign to flip on the log-determinant. That
is the whole difference from
[`mvg_pieces()`](https://statmodels7.github.io/distributions7/reference/mvg_pieces.md),
whose parametrization may carry either side.

## Notation

\\\mu\\ is the location, \\\Sigma\\ the scale matrix, \\\nu\\ the
degrees of freedom and \\\eta\\ the free vector of the matrix
parametrization.

## See also

[`mvstudent_t_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t_distrib.md)
for the family,
[`mvt_weights()`](https://statmodels7.github.io/distributions7/reference/mvt_weights.md)
for the quantities computed from these, and
[`mvg_pieces()`](https://statmodels7.github.io/distributions7/reference/mvg_pieces.md)
for the gaussian's twin.

## Examples

``` r
d <- mvstudent_t_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)
pc <- distributions7:::mvt_pieces(d, theta, derivs = TRUE)
names(pc)
#> [1] "mu"        "eta"       "nu"        "p"         "s"         "sigma"    
#> [7] "sigma_inv" "logdet"    "a"        

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
