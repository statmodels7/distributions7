# Multivariate Gaussian Generator

Draws from the multivariate gaussian by the standard factorization, \\y
= \mu + L z\\ with \\z\\ a vector of independent standard normals and
\\L L^\top = \Sigma\\. The factor is a Cholesky decomposition of the
covariance, taken after the precision has been inverted where the matrix
parametrization carries that side. The draws consume `n * p` values from
R's own generator through
[`stats::rnorm()`](https://rdrr.io/r/stats/Normal.html), so the stream
is reproducible under
[`base::set.seed()`](https://rdrr.io/r/base/Random.html).

## Arguments

- distrib:

  An
  [MvGaussianDistrib](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object, from
  [`mvgaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvgaussian1_distrib.md).

- n:

  The number of observations to draw. A single non-negative whole
  number; `n = 0` returns a matrix with zero rows.

- theta:

  A named list of parameters, each component a single number.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

An \\n \times p\\ numeric matrix, one draw per row, with column names
`v1`, ..., `vp` and no row names.

## See also

[`distrib_pdf.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.MvGaussianDistrib.md)
for the density this draws from,
[`mv_sigma.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MvGaussianDistrib.md)
for the covariance, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- mvgaussian1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4)

set.seed(1)
y <- distrib_rng(d, 4, theta)
y
#>              v1         v2
#> [1,] -0.1923385 -0.2808034
#> [2,]  0.7029573 -0.8982854
#> [3,] -0.4235124 -0.2351783
#> [4,]  2.2630579  0.9426015

# The sample moments approach the parameters, at the usual 1/sqrt(n) rate.
set.seed(2)
big <- distrib_rng(d, 20000, theta)
round(colMeans(big), 3)
#>     v1     v2 
#>  0.506 -0.299 
round(var(big), 3)
#>       v1    v2
#> v1 1.236 0.451
#> v2 0.451 0.839
round(mv_sigma(d, theta), 3)
#>       v1    v2
#> v1 1.221 0.442
#> v2 0.442 0.830
```
