# Marginal of a Multivariate Gaussian

Returns the marginal law of a subset of coordinates, which for a
gaussian is again a gaussian: the mean is the subvector and the
covariance the corresponding block, with nothing to integrate. The
result is a fresh
[`mvgaussian_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian_distrib.md)
of the reduced dimension, whose parameters are its own: a caller cannot
expect them to be a subset of the ones passed in.

## Arguments

- distrib:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian_distrib.md).

- theta:

  A named list of parameters, each component a single number.

- which:

  An integer vector of coordinates to keep, between 1 and \\p\\.
  Duplicates and out-of-range values are not checked and reach the
  matrix subsetting.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with `distrib`, an `MvGaussianDistrib` of dimension
`length(which)` on a log-Cholesky covariance, and `theta`, its
parameters as a named list.

## Details

The marginal is returned on an UNSTRUCTURED covariance whatever the
parent carried, because a block of a structured matrix need not have the
parent's structure. The leading block of an AR(1) is AR(1); a block of a
compound-symmetry matrix taken at scattered indices need not be; and a
block of a precision matrix is not the inverse of the corresponding
covariance block at all. Returning the unstructured form is correct in
every case, at the cost of `k(k+1)/2` free values where the parent may
have spent fewer.

## See also

[`mv_marginal.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.MvStudentTDistrib.md),
where the degrees of freedom are carried across unchanged,
[`plot.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/plot.multivariate_distrib.md),
whose panels are these marginals, and
[`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.md)
for the generic.

## Examples

``` r
d <- mvgaussian_distrib(3)
theta <- as.list(stats::setNames(
  c(1, -2, 0.5, 0.1, -0.2, 0.3, 0.4, -0.1, 0.2), d@params))

m <- mv_marginal(d, theta, c(1, 3))
m$distrib@n_dim
#> [1] 2

# The covariance block is carried across exactly.
all.equal(mv_sigma(m$distrib, m$theta),
          mv_sigma(d, theta)[c(1, 3), c(1, 3)], check.attributes = FALSE)
#> [1] TRUE

# And a single coordinate is the univariate normal, against stats::dnorm.
m1 <- mv_marginal(d, theta, 2)
s <- sqrt(mv_sigma(d, theta)[2, 2])
c(ours = distrib_pdf(m1$distrib, -1.4, m1$theta, log = TRUE),
  dnorm = dnorm(-1.4, mean = -2, sd = s, log = TRUE))
#>     ours    dnorm 
#> -1.04275 -1.04275 

# A precision parametrization marginalizes to a covariance, the block of
# the precision not being the precision of the block.
o <- mvgaussian_distrib(3, omega = parameters7::log_cholesky(3))
th_o <- as.list(stats::setNames(unlist(theta), o@params))
all.equal(mv_sigma(mv_marginal(o, th_o, c(1, 3))$distrib,
                   mv_marginal(o, th_o, c(1, 3))$theta),
          mv_sigma(o, th_o)[c(1, 3), c(1, 3)], check.attributes = FALSE)
#> [1] TRUE
```
