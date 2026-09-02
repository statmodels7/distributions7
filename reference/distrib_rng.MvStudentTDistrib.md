# Multivariate Student t Generator

Draws from the family through its scale-mixture representation, \\y =
\mu + L z\sqrt{\nu/g}\\ with \\z\\ a vector of independent standard
normals, \\g \sim \chi^2\_\nu\\ independent of them, and \\LL^\top =
\Sigma\\. A \\t\\ is a gaussian whose precision has been multiplied by a
gamma variate, and that same representation is what gives the family its
heavy tails and its closed-form expected information.

The draws consume `n * p` normal variates and `n` chi-squared variates
from R's own generators, in that order, so the stream is reproducible
under [`base::set.seed()`](https://rdrr.io/r/base/Random.html).

## Arguments

- distrib:

  An
  [MvStudentTDistrib](https://statmodels7.github.io/distributions7/reference/MvStudentTDistrib.md)
  object, from
  [`mvstudent_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/mvstudent_t1_distrib.md).

- n:

  The number of observations to draw. A single non-negative whole
  number.

- theta:

  A named list of parameters, each component a single number.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

An \\n \times p\\ numeric matrix, one draw per row, with column names
`v1`, ..., `vp`.

## Notation

\\\mu\\ is the location, \\\Sigma\\ the scale matrix, \\\nu\\ the
degrees of freedom and \\L\\ a lower Cholesky factor of \\\Sigma\\.

## See also

[`distrib_pdf.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.MvStudentTDistrib.md)
for the density this draws from,
[`variance.MvStudentTDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.MvStudentTDistrib.md)
for the covariance the sample approaches when it exists, and
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
for the generic.

## Examples

``` r
d <- mvstudent_t1_distrib(2)
theta <- list(mu1 = 0.5, mu2 = -0.3, sigma_log_L1 = 0.1,
              sigma_log_L2 = -0.2, sigma_L2.1 = 0.4, nu = 6)

set.seed(1)
distrib_rng(d, 4, theta)
#>              v1         v2
#> [1,] -0.1415980 -0.2822103
#> [2,]  0.7460944 -1.0254468
#> [3,] -0.8536934 -0.2049837
#> [4,]  2.9035115  1.3939925

# At nu = 6 the covariance exists and the sample approaches it.
set.seed(2)
big <- distrib_rng(d, 40000, theta)
round(var(big), 2)
#>      v1   v2
#> v1 1.85 0.67
#> v2 0.67 1.25
round(variance(d, theta), 2)
#>      v1   v2
#> v1 1.83 0.66
#> v2 0.66 1.25

# At nu = 1.5 it does not, and the sample covariance is a number that means
# nothing: it grows with the sample instead of settling.
t2 <- theta; t2$nu <- 1.5
set.seed(3)
vapply(c(2000, 20000, 200000),
       function(m) var(distrib_rng(d, m, t2))[1, 1], numeric(1))
#> [1]  38.15766 111.09339 406.30823
```
