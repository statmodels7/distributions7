# Derivatives of the Covariance with Respect to Every Parameter

Returns \\\partial\Sigma/\partial\theta_k\\ for EVERY parameter of a
multivariate distribution, in `distrib@params` order, so that the list
can be indexed by parameter position without any bookkeeping at the call
site. The parametrization supplies the derivatives of its own free
values; the location components and, for a Student \\t\\, the degrees of
freedom leave the matrix alone and contribute a matrix of zeros.

## Usage

``` r
mv_sigma_derivs(distrib, theta, n_before)
```

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object with a `param` property.

- theta:

  A named list of parameters, already aligned.

- n_before:

  The number of parameters before the matrix parametrization's free
  values, which is \\p\\ for both shipped families. Those positions get
  zero matrices.

## Value

A list of \\p \times p\\ numeric matrices, as long as `distrib@params`.
Every entry outside the parametrization's stretch is a matrix of zeros.

## Details

Where the parametrization carries the PRECISION the arrays are carried
onto the covariance by \\\partial\Sigma/\partial\eta_k = -\Sigma A_k
\Sigma\\, so the result is always in the covariance whichever side the
model is written on.

## Notation

\\\Sigma\\ is the covariance or scale matrix, \\\eta\\ the free vector
of the matrix parametrization, \\\theta\\ the distribution's parameter
vector and \\A_k\\ the parametrization's own derivative array.

## See also

[`mv_sd_cor()`](https://statmodels7.github.io/distributions7/reference/mv_sd_cor.md),
the consumer, and
[`mvg_pieces()`](https://statmodels7.github.io/distributions7/reference/mvg_pieces.md),
which does the same conversion for the gaussian's own methods.

## Examples

``` r
d <- mvgaussian_distrib(2)
theta <- distributions7:::align_theta(
  d, list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
          sigma_L2.1 = 0.5))
a <- distributions7:::mv_sigma_derivs(d, theta, 2)

# One entry per parameter, and the two mean entries are zero.
length(a) == d@n_params
#> [1] TRUE
a[[1]]
#> NULL

# The third entry is the derivative in the first free value, against a
# difference of mv_sigma().
h <- 1e-5
tp <- theta; tp$sigma_log_L1 <- tp$sigma_log_L1 + h
tm <- theta; tm$sigma_log_L1 <- tm$sigma_log_L1 - h
max(abs(a[[3]] - (mv_sigma(d, tp) - mv_sigma(d, tm)) / (2 * h)))
#> [1] 1.185736e-10
```
