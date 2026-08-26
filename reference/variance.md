# Variance of a Distribution or Sample

Computes \\\operatorname{Var}(Y) = E\[(Y - E\[Y\])^2\]\\ for a
distribution object, or the sample variance for a numeric vector.
Dispatch is on the first argument: a `distrib` uses the family's closed
form where one is registered and a quadrature or series otherwise, and a
numeric vector is passed to
[`stats::var()`](https://rdrr.io/r/stats/cor.html). Where the variance
does not exist the value is `NaN` or `Inf`, never a truncated
quadrature.

## Usage

``` r
variance(x, ...)
```

## Arguments

- x:

  An object inheriting from `distrib`, or a numeric vector.

- ...:

  For a `distrib`: `theta`, a named list of parameters, followed by any
  further arguments for
  [`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md).
  For a numeric vector: `na.rm`, a single logical, `FALSE` by default.

## Value

A numeric vector for a `distrib`, one value per parameter setting; a
single number for a numeric vector. `NaN` for a family with no moments,
`Inf` for one whose variance diverges.

## The two routes

\$\$\operatorname{Var}(Y) = \mathbb{E}\left\[(Y -
\mathbb{E}\[Y\])^{2}\right\].\$\$

43 of the 45 shipped families register a closed form, so the second
central moment is a formula in the parameters and costs about 24
microseconds; the two von Mises families fall through to
[`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md),
which is a quadrature and costs about 12 ms. The two routes agree to the
rounding of the quadrature, 3e-15 relative on a Gaussian.

## What the sample method returns

On a numeric vector the value is
[`stats::var()`](https://rdrr.io/r/stats/cor.html), which divides by
\\n - 1\\. That is the unbiased estimator, and it differs from the
convention
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
and
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)
use for their sample versions, which divide by \\n\\. The difference is
deliberate: both follow the commonest convention for the quantity in
question.

## See also

[`std_dev()`](https://statmodels7.github.io/distributions7/reference/std_dev.md)
for its square root,
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
and
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)
for the standardized third and fourth moments,
[`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md)
for the numerical route,
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
for the quadrature and the series.

## Examples

``` r
# A closed form, one value per setting.
variance(gaussian1_distrib(), list(mu = 0, sigma = c(1, 2, 4)))
#> [1]  1  4 16

# A Poisson is equidispersed: the variance is the mean.
variance(poisson_distrib(), list(mu = 3))
#> [1] 3

# A Student t has a variance only above two degrees of freedom.
variance(student_t1_distrib(), list(mu = 0, sigma = 1, nu = c(3, 10)))
#> [1] 3.00 1.25

# On a numeric vector this is stats::var, with the n - 1 denominator.
set.seed(1)
y <- rnorm(50)
all.equal(variance(y), var(y))
#> [1] TRUE
```
