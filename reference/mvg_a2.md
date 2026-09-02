# The Second Derivative of the Covariance, by Position

Reads the component of
[`parameters7::param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.html)
belonging to the pair of free values \\(k, l)\\. The key is constructed
from the sorted pair of free names, never parsed out of a stored key, so
a free value whose own label contains the separator cannot split into
the wrong number of pieces. The result is symmetric in its two
positions.

## Usage

``` r
mvg_a2(pc, k, l)
```

## Arguments

- pc:

  The result of
  [`mvg_pieces()`](https://statmodels7.github.io/distributions7/reference/mvg_pieces.md)
  called with `derivs2 = TRUE`, whose `a2` component holds the second
  derivative arrays and whose `s` component supplies the free names.

- k, l:

  Positions among the matrix parametrization's free values, each a
  single whole number between 1 and `s@n_free`. Their order does not
  matter.

## Value

A \\p \times p\\ symmetric numeric matrix,
\\\partial^2\Sigma/\partial\eta_k\partial\eta_l\\.

## See also

[`mvg_pieces()`](https://statmodels7.github.io/distributions7/reference/mvg_pieces.md)
for the argument and
[`param_pair_lookup()`](https://statmodels7.github.io/distributions7/reference/param_pair_lookup.md)
for the other route to the same enumeration.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)
pc <- distributions7:::mvg_pieces(d, theta, derivs2 = TRUE)

distributions7:::mvg_a2(pc, 1L, 3L)
#>          [,1]     [,2]
#> [1,] 0.000000 1.105171
#> [2,] 1.105171 0.000000

# Symmetric in the two positions.
identical(distributions7:::mvg_a2(pc, 1L, 3L),
          distributions7:::mvg_a2(pc, 3L, 1L))
#> [1] TRUE

# And it is the mixed second derivative, against a difference of the first.
h <- 1e-5
eta <- c(0.1, -0.2, 0.4)
s <- d@param
num <- (parameters7::param_d1(s, eta + c(0, 0, h))[[1]] -
        parameters7::param_d1(s, eta - c(0, 0, h))[[1]]) / (2 * h)
max(abs(distributions7:::mvg_a2(pc, 1L, 3L) - num))
#> [1] 2.655653e-13
```
