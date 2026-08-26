# Skewness of a Distribution or Sample

Computes the skewness, the third standardized central moment

\$\$\gamma_1 = \mathbb{E}\\\left\[\left(\frac{Y -
\mathbb{E}\[Y\]}{\operatorname{sd}(Y)}\right)^{3}\right\],\$\$

for a distribution object, or the sample skewness for a numeric vector.
It is zero for a symmetric law, positive for a right tail and negative
for a left one. 43 of the 45 shipped families register a closed form;
the two von Mises reach a quadrature through
[`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md).

## Usage

``` r
skewness(x, ...)
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
`Inf` where the third moment diverges.

## What the sign and the size mean

\\\gamma_1\\ is invariant to location and to positive scaling, so it is
a property of the shape alone: every Gaussian has \\\gamma_1 = 0\\ and
every exponential has \\\gamma_1 = 2\\, whatever their parameters. A
family with a shape parameter moves along the axis with it, and that is
why the quantity is worth reporting: a Poisson at mean \\\mu\\ has
\\\gamma_1 = \mu^{-1/2}\\, so its asymmetry vanishes as the counts grow.

## When it fails to exist

The third moment is needed, so a family with heavy tails may have no
skewness even where it has a mean. A Student t has one only above three
degrees of freedom and returns `Inf` at or below them; a Cauchy has none
at any parameter value and returns `NaN`.

## What the sample method returns

On a numeric vector the third and second central moments both divide by
\\n\\, so the value is \\m_3 / m_2^{3/2}\\ with \\m_k = n^{-1}\sum_i
(y_i - \bar y)^k\\. This is the population denominator, and it differs
from the \\n - 1\\ convention
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
uses on a vector.

## See also

[`kurtosis()`](https://statmodels7.github.io/distributions7/reference/kurtosis.md)
for the fourth standardized moment,
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
and
[`std_dev()`](https://statmodels7.github.io/distributions7/reference/std_dev.md)
for the second,
[`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md)
for the numerical route.

## Examples

``` r
# Zero for a symmetric family, at every parameter value.
skewness(gaussian1_distrib(), list(mu = c(-2, 0, 5), sigma = c(1, 2, 3)))
#> [1] 0 0 0

# A shape parameter moves it: a Poisson is mu^(-1/2).
all.equal(skewness(poisson_distrib(), list(mu = 4)), 0.5)
#> [1] TRUE

# A gamma is 2 / sqrt(shape), so it flattens as the shape grows.
skewness(gamma2_distrib(), list(mu = 2, sigma2 = c(4, 1, 0.25)))
#> [1] 2.0 1.0 0.5

# The sample version uses the n denominator.
set.seed(1)
y <- rgamma(200, shape = 2)
c(sample = skewness(y), theory = 2 / sqrt(2))
#>    sample    theory 
#> 0.8573301 1.4142136 
```
