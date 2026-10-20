# A Proposal for Integrating a Multivariate Density

Draws from a distribution that dominates the family, and returns the
log-density of the draws under that distribution, so that the
importance-sampling estimate of \\\int f = 1\\ can be formed.
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
is its only consumer; a discrete family takes the exact sum over
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

  A named list or vector of parameters, in which each component is a
  single number.

- n:

  The number of draws, a single positive whole number.

- ...:

  Passed to methods. The shipped methods do not read it.

## Value

A named list with `y`, an \\n \times p\\ numeric matrix of draws, and
`logd`, a numeric vector of length \\n\\ holding their log-density under
the proposal.

## The default proposal

The method of the class draws from a gaussian with the family's mean and
twice its covariance, which requires a finite, non-singular covariance.
The inflation keeps the estimate informative: for a gaussian family, a
proposal equal to the family makes every importance ratio one.

## Families with a proposal of their own

A Dirichlet lives on the simplex, a set of measure zero in
\\\mathbb{R}^p\\, where a gaussian spread over the ambient space places
no mass; its covariance is singular, and the estimate from the default
proposal is far from one. The Dirichlet registers the uniform
distribution on the simplex, whose density is the constant
\\\Gamma(p)\\. A Student t has no finite covariance at \\\nu \le 2\\,
and its tails are heavier than a gaussian's at every \\\nu\\; it
registers a Student t with half its degrees of freedom and twice its
scale matrix.

The draws and their log-density are taken against the dominating measure
with respect to which the family's density is written.

## Notation

\\f\\ is the family's density, \\p\\ the dimension and \\\Gamma\\ the
gamma function.

## See also

[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md),
the only consumer,
[`mv_reference_draw.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/mv_reference_draw.multivariate_distrib.md)
for the default,
[`mv_reference_draw.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_reference_draw.MvStudentTDistrib.md)
for the Student t's, and
[`mv_support()`](https://statmodels7.github.io/distributions7/reference/mv_support.md),
the exact route of a discrete family.

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

# This is the estimate for which the proposal is built: the density over
# the proposal, averaged.
set.seed(2)
r <- mv_reference_draw(d, theta, 20000)
mean(exp(distrib_pdf(d, r$y, theta, log = TRUE) - r$logd))
#> [1] 0.992289

# A Dirichlet registers its own proposal, the uniform distribution on the
# simplex.
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
