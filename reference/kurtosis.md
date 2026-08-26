# Excess Kurtosis of a Distribution or Sample

Computes the excess kurtosis, the fourth standardized central moment
less the three a Gaussian has:

\$\$\gamma_2 = \mathbb{E}\\\left\[\left(\frac{Y -
\mathbb{E}\[Y\]}{\operatorname{sd}(Y)}\right)^{4}\right\] - 3.\$\$

The subtraction puts the Gaussian at zero, so the sign reads as a
comparison with it: positive for heavier tails and a sharper peak,
negative for lighter ones. 42 of the 45 shipped families register a
closed form; the elastic net and the two von Mises reach a quadrature
through
[`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md).

## Usage

``` r
kurtosis(x, ...)
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
`Inf` where the fourth moment diverges.

## Why the three is subtracted

The raw fourth standardized moment is 3 for every Gaussian, whatever its
mean and scale, so it carries no information about the Gaussian itself.
Subtracting it makes the quantity a signed distance from normality and
puts the two commonest reference laws at recognizable values: 0 for the
Gaussian, 3 for the Laplace, \\6/(\nu-4)\\ for a Student t.

## When it fails to exist

The fourth moment is needed, so the threshold is higher than for
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md):
a Student t has an excess kurtosis only above four degrees of freedom
and returns `Inf` at or below them, and a Cauchy returns `NaN` at every
parameter value.

## What the sample method returns

On a numeric vector the fourth and second central moments both divide by
\\n\\, giving \\m_4 / m_2^{2} - 3\\. This is the population denominator,
the same convention
[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
uses.
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
divides by \\n - 1\\ instead.

## See also

[`skewness()`](https://statmodels7.github.io/distributions7/reference/skewness.md)
for the third standardized moment,
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
and
[`std_dev()`](https://statmodels7.github.io/distributions7/reference/std_dev.md)
for the second,
[`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md)
for the numerical route.

## Examples

``` r
# Zero for every Gaussian, by construction.
kurtosis(gaussian1_distrib(), list(mu = c(-2, 0, 5), sigma = c(1, 2, 3)))
#> [1] 0 0 0

# Three for the Laplace, which is the textbook heavy-tailed comparison.
all.equal(kurtosis(laplace_distrib(), list(mu = 0, sigma = 2)), 3)
#> [1] TRUE

# A Student t is 6 / (nu - 4), and infinite at or below four.
kurtosis(student_t1_distrib(), list(mu = 0, sigma = 1, nu = c(3, 5, 10)))
#> [1] Inf   6   1

# The sample version uses the n denominator.
set.seed(1)
y <- rt(500, df = 10)
c(sample = kurtosis(y), theory = 6 / (10 - 4))
#>   sample   theory 
#> 1.271059 1.000000 
```
