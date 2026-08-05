# A Proposal for Integrating a Multivariate Density

Draws from a distribution that dominates the family, together with the
log-density of those draws, so that an importance-sampling estimate of
\\\int f = 1\\ can be formed.

## Usage

``` r
mv_reference_draw(distrib, theta, n, ...)
```

## Arguments

- distrib:

  An object inheriting from class
  [`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md).

- theta:

  A named list or vector of parameters.

- n:

  The number of draws.

- ...:

  Passed to methods.

## Value

A list with `y`, a matrix of `n` draws, and `logd`, their log-density
under the proposal.

## Details

The default proposal is a gaussian with the same mean and twice the
covariance, which serves any family supported on all of
\\\mathbb{R}^p\\. The inflation matters: a proposal equal to the
distribution itself makes every ratio one and certifies nothing.

A family whose support is a lower-dimensional subset of
\\\mathbb{R}^p\\, such as a Dirichlet on the simplex, must register its
own method, because a proposal spread over the ambient space places no
mass at all on the support and the estimate would be zero. The draws and
the log-density must be taken with respect to the same dominating
measure the family's density is written against.

This is consumed by
[`check_distrib`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
and by nothing else. A discrete family does not need it, its
normalization being an exact sum over
[`mv_support`](https://statmodels7.github.io/distributions7/reference/mv_support.md).

## See also

[`check_distrib`](https://statmodels7.github.io/distributions7/reference/check_distrib.md),
[`mv_support`](https://statmodels7.github.io/distributions7/reference/mv_support.md)

## Examples

``` r
d <- mvgaussian_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0)
set.seed(1)
str(mv_reference_draw(d, theta, 5))
#> List of 2
#>  $ y   : num [1:5, 1:2] -0.886 0.26 -1.182 2.256 0.466 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:2] "v1" "v2"
#>  $ logd: num [1:5] -3.06 -2.67 -3.15 -3.97 -2.63
```
