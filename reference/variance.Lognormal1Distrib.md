# Variance of the Lognormal Distribution

Closed form: \$\$\operatorname{Var}(Y) = (e^{\sigma^2} - 1)\\ e^{2\mu +
\sigma^2}.\$\$ The two factors separate the two effects:
\\e^{2\mu+\sigma^2}\\ is the square of the mean and sets the scale, and
\\e^{\sigma^2} - 1\\ is the squared coefficient of variation and depends
on the spread of the logarithm alone.

## Arguments

- x:

  A `Lognormal1Distrib`, from
  [`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md).

- theta:

  A named list with components `mu` (any real value) and `sigma2`
  (positive), each a numeric vector of length 1 or `n`. The variance
  grows like \\e^{2\sigma^2}\\ and overflows to `Inf` well before
  `sigma2` reaches 400.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of variances, of length
`max(length(theta$mu), length(theta$sigma2))`.

## Notation

\\\mu\\ is the mean of \\\log Y\\ and \\\sigma^2 \> 0\\ its variance.

## See also

[`mean.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Lognormal1Distrib.md),
[`skewness.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Lognormal1Distrib.md),
[`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md).

## Examples

``` r
d <- lognormal1_distrib()

# (exp(sigma2) - 1) exp(2 mu + sigma2).
all.equal(variance(d, list(mu = 0, sigma2 = 1)), (exp(1) - 1) * exp(1))
#> [1] TRUE

# The coefficient of variation depends on sigma2 alone.
th <- list(mu = c(-2, 0, 5), sigma2 = 0.25)
std_dev(d, th) / mean(d, th)
#> [1] 0.5329404 0.5329404 0.5329404
```
