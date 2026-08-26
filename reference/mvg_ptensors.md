# Precision Derivative Arrays of a Multivariate Gaussian

Supplies the derivative arrays of the PRECISION \\P = \Sigma^{-1}\\ in
the matrix parametrization's free values, to orders 1 through 4, as an
accessor keyed by an index multiset. Under a precision parametrization
these are the parametrization's own `param_d1()` to `param_d4()`; under
a covariance one they follow from repeated differentiation of the
inverse, \$\$P_t = \sum\_{(B_1,\dots,B_q)} (-1)^q\\ P A\_{B_1} P \cdots
A\_{B_q} P,\$\$ summed over the ordered partitions of the multiset \\t\\
into nonempty blocks, with \\A_B\\ the derivative of \\\Sigma\\ in the
free values \\B\\. Nothing is transcribed from an expanded formula, so
no term can go missing at order 3 or 4.

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
  array up to that order is enumerated.

- inverted:

  Logical of length 1. `TRUE` when the free values parametrize the
  precision, in which case the arrays are read off the parametrization
  directly and no expansion runs. Defaults to `FALSE`.

## Value

A named list with `get`, a function of one integer vector returning the
\\p \times p\\ array for that multiset; `sign_ld`, the coefficient the
log-determinant term carries in a derivative of the log-density
(\\-1/2\\ for a covariance, \\+1/2\\ for a precision); and `pc`, the
pieces as supplied.

## Details

The accessor memoizes, so an array asked for twice within one call is
computed once. It answers for the EMPTY multiset with \\P\\ itself: the
gaussian never asks for that, but the multivariate Student t does, a
partition block of pure mean indices carrying no matrix index at all.

The function takes the pieces, so that the multivariate Student t can
hand it the same arrays of its scale matrix and the toolkit carries one
copy of this expansion. That copy's first draft double counted the mixed
terms, and only a comparison against a finite-difference stencil caught
it.

## Notation

\\\Sigma\\ is the covariance, \\P = \Sigma^{-1}\\ the precision,
\\\eta\\ the free vector of the matrix parametrization, \\A_t\\ its
derivative array for the multiset \\t\\, and \\P_t\\ the corresponding
derivative of the precision.

## See also

[`mv_ordered_partitions()`](https://statmodels7.github.io/distributions7/reference/mv_ordered_partitions.md)
for the enumeration it sums over,
[`mvg_higher()`](https://statmodels7.github.io/distributions7/reference/mvg_higher.md)
for the consumer, and
[`mvg_pieces()`](https://statmodels7.github.io/distributions7/reference/mvg_pieces.md)
for the argument.

## Examples

``` r
d <- mvgaussian_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
pc <- distributions7:::mvg_pieces(d, theta, derivs2 = TRUE)
pt <- distributions7:::mvg_ptensors(pc, 2L)

P <- pc$sigma_inv
A1 <- pc$a[[1]]

# The empty multiset is the precision itself.
all.equal(pt$get(integer(0)), P)
#> [1] TRUE

# At order one the expansion is the derivative of an inverse.
all.equal(pt$get(1L), -P %*% A1 %*% P)
#> [1] TRUE

# At order two it is the three-term expansion, written out here.
A3 <- pc$a[[3]]
A13 <- distributions7:::mvg_a2(pc, 1L, 3L)
all.equal(pt$get(c(1L, 3L)),
          P %*% A1 %*% P %*% A3 %*% P + P %*% A3 %*% P %*% A1 %*% P -
            P %*% A13 %*% P)
#> [1] TRUE

# A precision parametrization needs no expansion and flips the sign the
# log-determinant term carries.
c(covariance = pt$sign_ld,
  precision = distributions7:::mvg_ptensors(pc, 2L, inverted = TRUE)$sign_ld)
#> covariance  precision 
#>       -0.5        0.5 
```
