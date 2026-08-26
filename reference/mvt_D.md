# The Log-Ratio Cancellation of a Multivariate t

Computes \\D(u) = \tfrac{u}{1+u} - \log(1+u)\\, the data-carrying part
of the score in \\\nu\\, taken at \\u = q/\nu\\. The two terms agree to
first order, so below \\\lvert u\rvert = 10^{-3}\\ the expansion
\\-u^2/2 + 2u^3/3 - 3u^4/4 + 4u^5/5\\ is used instead of the difference.

## Usage

``` r
mvt_D(u)
```

## Arguments

- u:

  A numeric vector of any length, taken as \\q/\nu\\. Values at or below
  \\-1\\ are outside the domain of the logarithm and give `NaN` from
  [`base::log1p()`](https://rdrr.io/r/base/Log.html); the intended
  argument is non-negative.

## Value

A numeric vector as long as `u`, non-positive, of order \\-u^2/2\\ near
zero.

## Details

The argument \\u\\ is the squared Mahalanobis distance divided by the
degrees of freedom, so it is small exactly where the family is close to
the gaussian or the observation is close to the location, so it is the
ordinary case and not an edge one. Taken directly at \\u = 10^{-8}\\ the
difference loses about nine significant figures.

## Notation

\\q\\ is the squared Mahalanobis distance of an observation from the
location and \\\nu\\ the degrees of freedom.

## See also

[`mvt_A()`](https://statmodels7.github.io/distributions7/reference/mvt_A.md)
and
[`mvt_T()`](https://statmodels7.github.io/distributions7/reference/mvt_T.md)
for the other two cancellations, and
[`distrib_gradient.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MvStudentTDistrib.md)
for the consumer.

## Examples

``` r
u <- c(1e-8, 1e-5, 1e-3, 0.1, 1)

# Near zero the value is -u^2 / 2, and it approaches that ratio smoothly.
cbind(u = u, scaled = distributions7:::mvt_D(u) / (-u^2 / 2))
#>          u    scaled
#> [1,] 1e-08 1.0000000
#> [2,] 1e-05 0.9999867
#> [3,] 1e-03 0.9986682
#> [4,] 1e-01 0.8802178
#> [5,] 1e+00 0.3862944

# Against the direct difference, relatively. The two agree where the
# difference still has digits and part company as u shrinks.
cbind(u = u,
      rel_gap = abs(distributions7:::mvt_D(u) - (u / (1 + u) - log1p(u))) /
        abs(distributions7:::mvt_D(u)))
#>          u      rel_gap
#> [1,] 1e-08 2.891508e-09
#> [2,] 1e-05 6.023506e-12
#> [3,] 1e-03 0.000000e+00
#> [4,] 1e-01 0.000000e+00
#> [5,] 1e+00 0.000000e+00
```
