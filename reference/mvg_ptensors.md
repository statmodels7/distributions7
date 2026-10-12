# Precision Derivative Arrays of a Multivariate Gaussian

Supplies the derivative arrays of the precision \\P = \Sigma^{-1}\\ in
the free values of the matrix parametrization, to orders 1 through 4, as
an accessor keyed by an index multiset. Under a precision
parametrization they are the parametrization's own `param_d1()` to
`param_d4()`; under a covariance parametrization they are its
derivatives of the inverse,
[`parameters7::param_inv_d1()`](https://statmodels7.github.io/parameters7/reference/param_inv_d1.html)
to
[`parameters7::param_inv_d4()`](https://statmodels7.github.io/parameters7/reference/param_inv_d1.html).
For
[`parameters7::log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.html)
these are exact where the covariance is nearly singular; for another
parametrization they are the ordered-partition sum of
[`parameters7::inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.html).

## Usage

``` r
mvg_ptensors(pc, order, inverted = FALSE)
```

## Arguments

- pc:

  The pieces, as returned by
  [`mvg_pieces()`](https://statmodels7.github.io/distributions7/reference/mvg_pieces.md)
  or
  [`mvt_pieces()`](https://statmodels7.github.io/distributions7/reference/mvt_pieces.md):
  any list carrying `s` (the parametrization), `eta` (its free vector)
  and `sigma_inv`.

- order:

  The highest order wanted, a single whole number from 1 to 4. Every
  array up to that order is computed.

- inverted:

  Logical of length 1. `TRUE` when the free values parametrize the
  precision. Defaults to `FALSE`.

## Value

A named list with `get`, a function of one integer vector returning the
\\p \times p\\ array for that multiset; `sign_ld`, the coefficient that
the log-determinant term carries in a derivative of the log-density
(\\-1/2\\ for a covariance, \\+1/2\\ for a precision); and `pc`, the
pieces as supplied.

## Details

The accessor returns \\P\\ itself for the empty multiset. The gaussian
does not use it; the multivariate Student t does, for a partition block
of mean indices that carries no matrix index. The function takes the
pieces, so that the Student t passes it the arrays of its scale matrix
and the package keeps one copy of the construction.

## Notation

\\\Sigma\\ is the covariance, \\P = \Sigma^{-1}\\ the precision,
\\\eta\\ the free vector of the matrix parametrization and \\P_t\\ the
derivative of the precision over the multiset \\t\\.

## See also

[`mvg_higher()`](https://statmodels7.github.io/distributions7/reference/mvg_higher.md)
and
[`mvt_higher()`](https://statmodels7.github.io/distributions7/reference/mvt_higher.md)
for the consumers, and
[`mvg_pieces()`](https://statmodels7.github.io/distributions7/reference/mvg_pieces.md)
for the argument.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
pc <- distributions7:::mvg_pieces(d, theta, derivs2 = TRUE)
pt <- distributions7:::mvg_ptensors(pc, 2L)

P <- pc$sigma_inv
A1 <- pc$a[[1]]

# The empty multiset is the precision itself.
all.equal(pt$get(integer(0)), P)
#> [1] TRUE

# At order one it is the derivative of an inverse, -P A P.
all.equal(pt$get(1L), -P %*% A1 %*% P)
#> [1] TRUE

# A precision parametrization reads its own arrays and flips the sign that
# the log-determinant term carries.
c(covariance = pt$sign_ld,
  precision = distributions7:::mvg_ptensors(pc, 2L, inverted = TRUE)$sign_ld)
#> covariance  precision 
#>       -0.5        0.5 
```
