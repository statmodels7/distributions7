# Standard Deviation of a Distribution or Sample

Computes \\\operatorname{sd}(Y) = \sqrt{\operatorname{Var}(Y)}\\ for a
distribution object, or the sample standard deviation for a numeric
vector. No family registers a closed form of its own, so a `distrib`
always reaches the square root of
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md),
which is where a family's formula is consulted. A numeric vector is
passed to [`stats::sd()`](https://rdrr.io/r/stats/sd.html).

## Usage

``` r
std_dev(x, ...)
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
single number for a numeric vector.

## Details

\$\$\operatorname{sd}(Y) = \sqrt{\operatorname{Var}(Y)}.\$\$

The square root is taken after the variance, so accuracy and cost are
the variance's: a closed form for 43 of the 45 shipped families, a
quadrature for the two von Mises. A family whose variance is `NaN` or
`Inf` gives the same here, and a negative variance cannot arise, so the
root is always real.

On a numeric vector the value is
[`stats::sd()`](https://rdrr.io/r/stats/sd.html), the root of the \\n -
1\\ variance.

## See also

[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md),
of which this is the square root;
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
and
[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md),
which standardize by it;
[`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md)
for the numerical route.

## Examples

``` r
# The scale of a Gaussian is its standard deviation, exactly.
std_dev(gaussian1_distrib(), list(mu = 0, sigma = c(1, 2, 4)))
#> [1] 1 2 4

# The square root of the variance, on any family.
d <- gamma2_distrib()
all.equal(std_dev(d, list(mu = 2, sigma2 = 1)),
          sqrt(variance(d, list(mu = 2, sigma2 = 1))))
#> [1] TRUE

# On a numeric vector this is stats::sd.
set.seed(1)
y <- rnorm(50)
all.equal(std_dev(y), sd(y))
#> [1] TRUE
```
