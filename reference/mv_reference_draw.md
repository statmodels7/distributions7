# A Proposal for Integrating a Multivariate Density

Draws from a distribution that DOMINATES the family, together with the
log-density of those draws, so that an importance-sampling estimate of
\\\int f = 1\\ can be formed. It is consumed by
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
and by nothing else, a discrete family taking the exact sum over
[`mv_support()`](https://statmodels7.github.io/distributions7/reference/mv_support.md)
instead.

## Usage

``` r
mv_reference_draw(distrib, theta, n, ...)
```

## Arguments

- distrib:

  An object inheriting from
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md).

- theta:

  A named list or vector of parameters, each component a single number.

- n:

  The number of draws, a single positive whole number.

- ...:

  Passed to methods. No shipped method reads it.

## Value

A named list with `y`, an \\n \times p\\ numeric matrix of draws, and
`logd`, a numeric vector of length \\n\\ holding their log-density under
the PROPOSAL.

## Why the proposal is inflated

The default is a gaussian with the family's mean and TWICE its
covariance, which serves any family supported on all of
\\\mathbb{R}^p\\. The inflation is load-bearing: a proposal equal to the
distribution itself makes every importance ratio one, so the estimate is
one by construction and certifies nothing.

## A family on a lower-dimensional support must override

A Dirichlet lives on the simplex, a set of measure zero in
\\\mathbb{R}^p\\, and a gaussian proposal spread over the ambient space
places no mass on it. The failure is QUIET:
[`chol()`](https://rdrr.io/r/base/chol.html) accepts the singular
covariance, and the estimate of an integral that is 1 comes back at
\\2\times10^{-8}\\. The Dirichlet therefore registers the uniform on the
simplex, whose density is the constant \\\Gamma(p)\\.

The draws and the log-density must be taken against the SAME dominating
measure the family's density is written against, or the ratio means
nothing.

## Notation

\\f\\ is the family's density, \\p\\ the dimension and \\\Gamma\\ the
gamma function.

## See also

[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md),
the only consumer,
[`mv_reference_draw.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/mv_reference_draw.multivariate_distrib.md)
for the default, and
[`mv_support()`](https://statmodels7.github.io/distributions7/reference/mv_support.md),
the exact route a discrete family takes.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0, mu2 = 0, sigma_log_L1 = 0, sigma_log_L2 = 0,
              sigma_L2.1 = 0)
set.seed(1)
str(mv_reference_draw(d, theta, 5))
#> List of 2
#>  $ y   : num [1:5, 1:2] -0.886 0.26 -1.182 2.256 0.466 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:2] "v1" "v2"
#>  $ logd: num [1:5] -3.06 -2.67 -3.15 -3.97 -2.63

# The estimate it is built for: the density over the proposal, averaged.
set.seed(2)
r <- mv_reference_draw(d, theta, 20000)
mean(exp(distrib_pdf(d, r$y, theta, log = TRUE) - r$logd))
#> [1] 0.992289

# A Dirichlet lives on the simplex and registers its own proposal, the
# uniform there.
dd <- dirichlet_distrib(3)
thd <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)
set.seed(3)
rd <- mv_reference_draw(dd, thd, 5)
round(rd$y, 4)
#>        [,1]   [,2]   [,3]
#> [1,] 0.8582 0.1035 0.0383
#> [2,] 0.1860 0.6907 0.1232
#> [3,] 0.8799 0.0072 0.1129
#> [4,] 0.1896 0.0129 0.7975
#> [5,] 0.2266 0.1271 0.6463
unique(round(rd$logd, 10))
#> [1] 0.6931472
log(gamma(3))
#> [1] 0.6931472
```
